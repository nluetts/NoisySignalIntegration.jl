const PRIMARY_COLOR = :blue
const SECONDARY_COLOR = :red
OFFSET_FIT_PLOT = 0.75
WIDTH_FIT_PLOT = 4.0


# this function is required to plot integration areas
function get_left_right_points(
    xs::AbstractArray{T},
    ys::AbstractArray{T},
    xₗ::T,
    xᵣ::T,
    b::Union{Nothing,UncertainBound}
    ;
    baseline_handling=nothing
) where {T<:AbstractFloat}

    xₗ, xᵣ = xₗ < xᵣ ? (xₗ, xᵣ) : (xᵣ, xₗ)

    i = searchsortedfirst(xs, xₗ)
    j = searchsortedfirst(xs, xᵣ)
    any(
        [
        i < 2,
        i > length(xs) - 1,
        j < 3,
        j > length(xs)
    ]
    ) && throw(error("At least one integration bound is outside the support region ($(minimum(xs)), $(maximum(xs)))."))


    if baseline_handling == "local"
        _, _, yₗ, yᵣ = _local_baseline(xs, ys, xₗ, xᵣ, b)
    else
        # bound `b` is unused in this case
        _, _, yₗ, yᵣ = _endpoint_to_endpoint_baseline(xs, ys, xₗ, xᵣ)
    end

    return i, j, xₗ, xᵣ, yₗ, yᵣ
end

@recipe function plot_recipe(crv::AbstractCurve)
    return crv.x, crv.y
end


# --------------------------------------
# enable plotting of sample draws
# --------------------------------------

@recipe function plot_recipe(c::AbstractCurve, uc::UncertainCurve, draws=3)
    draws < 0 && throw(ArgumentError("Number of samples must be > 0."))

    legend := :none
    layout := (draws + 1, 1)
    link := :both
    size --> (500, 600)

    delete!(plotattributes, :draws)

    for i ∈ 0:draws
        @series begin
            subplot := i + 1
            if i == 0
                seriescolor := SECONDARY_COLOR
                yguide := "input"
                c
            else
                yguide := "sample $(i)"
                get_draw(i, uc)
            end
        end
    end
end


@recipe function plot_recipe(crv::Curve{T},
    left::T,
    right::T; subtract_baseline=false,
    local_baseline=false,
    bound=nothing,
    draw_band_centers=false,
) where T

    (local_baseline && bound == nothing) && error("You have to provide a bound if local_baseline == true.") |> throw
    (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

    left = T(left)
    right = T(right)

    # draw area
    if subtract_baseline
        bh = "end-to-end"
    elseif local_baseline
        bh = "local"
    else
        bh = nothing
    end
    l, r, xl, xr, yl, yr = get_left_right_points(crv.x, crv.y, left, right, bound; baseline_handling=bh)

    if !(local_baseline || subtract_baseline)
        x = [xl; crv.x[l+1:r-1]; xr; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; zero(T); zero(T)]
    elseif subtract_baseline
        x = [xl; crv.x[l+1:r-1]; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; yl]
    else
        _, _, ycl, ycr = _endpoint_to_endpoint_baseline(crv.x, crv.y, left, right) # to get extra points located on the curve
        x = [xl; crv.x[l+1:r-1]; xr; xr; xl]
        y = [ycl; crv.y[l+1:r-1]; ycr; yr; yl]
    end

    @series begin
        fillrange := 0
        fillalpha --> 0.5
        fillcolor --> :orange
        linewidth --> 0
        label --> nothing
        x, y
    end

    if draw_band_centers
        # draw band center
        if subtract_baseline
            bc = band_center(crv, xl, xr, true)
            _, _, yl, yr = _endpoint_to_endpoint_baseline(crv.x, crv.y, xl, xr)
            y0 = lininterp(bc, xl, xr, yl, yr)
        elseif local_baseline
            baseline = baseline_from_points(_local_baseline(crv.x, crv.y, xl, xr, bound)...)
            bc = band_center(crv, xl, xr, baseline)
            _, _, yl, yr = _local_baseline(crv.x, crv.y, xl, xr, bound)
            y0 = lininterp(bc, xl, xr, yl, yr)
        else
            bc = band_center(crv, xl, xr)
            y0 = zero(typeof(crv.y[1]))
        end
        y = lininterp(pmean(bc), crv)
        @series begin
            alpha --> 0.3
            color --> :black
            label --> nothing
            [bc, bc], [y0, y]
        end
    end
end


@recipe function plot_recipe(crv::Curve{T}, bnds::Vector{UncertainBound{T,N}}, draw::Int) where {T,N}
    @series begin
        crv
    end
    for (j, b) ∈ enumerate(bnds)
        @series begin
            fillcolor := j % 2 == 1 ? :red : :orange
            left, right = get_draw(draw, b)
            bound := b
            crv, left, right
        end
    end
end


@recipe function plot_recipe(crv::Curve{T}, bnd::UncertainBound{T,N}, draw::Int) where {T,N}
    @series begin
        crv, [bnd], draw
    end
end


# plot draws of curves alongside with draws of bounds
@recipe function plot_recipe(
    uc::UncertainCurve{T,N},
    bnds::Vector{UncertainBound{T,N}}
    ;
    draws=3,
    subtract_baseline=true
) where {T,N}

    # Be flexible here, if user gives a range, we draw these specific draws
    padleft(range) = (max(0, minimum(draws)-1)):maximum(draws)
    range = draws isa UnitRange ? padleft(draws) : 0:draws

    legend := :none
    layout := (length(range), 1)
    link := :both
    size --> (500, 600)

    mean_uc = mean(uc)
    
    for (k, i) ∈ enumerate(range)
        for (j, b) in enumerate(bnds)
            @series begin
                fillcolor := j % 2 == 1 ? :red : :orange
                subplot := k
                bound := b
                if k == 1
                    mean_uc, mean(b)...
                else
                    get_draw(i, uc), get_draw(i, b)...
                end
            end
        end
        @series begin
            subplot := k
            if k == 1
                # mean spectrum
                seriescolor := SECONDARY_COLOR
                yguide := "mean"
                mean_uc
            else
                seriescolor := PRIMARY_COLOR
                yguide := "sample $(i)"
                get_draw(i, uc)
            end
        end
    end
end


@recipe function plot_recipe(
    uc::UncertainCurve{T,N},
    bnd::UncertainBound{T,N}
    ;
    draws=3
) where {T,N}
    draws := draws
    uc, [bnd]
end



# --------------------------------------
# enable plotting of noise sample draws
# --------------------------------------

@recipe function plot_recipe(x::Vector{T}, nm::AbstractNoiseModel; draws=3, subplot_offset=0) where {T<:Real}
    draws < 0 && throw(ArgumentError("Number of samples must be > 0."))

    layout --> (draws, 1)
    legend --> :none
    link --> :both
    size --> (500, 600)

    delete!(plotattributes, :grid_points)
    delete!(plotattributes, :draws)

    S = generate_noise(nm, length(x), draws)
    for i ∈ 1:draws
        @series begin
            subplot := i + subplot_offset # the offset does only apply if plotting together with a noise sample
            yguide := "sample $(i)"
            x, get_draw.(i, S)
        end
    end
end
@recipe function plot_recipe(nm::AbstractNoiseModel; gridpoints=1000, draws=3)
    @series begin
        x = collect(eltype(nm), 1:gridpoints)
        x, nm
    end
end


@recipe function plot_recipe(ns::NoiseSample, nm::AbstractNoiseModel; draws=3)

    layout --> (draws + 1, 1)
    legend --> :none
    link --> :both
    size --> (500, 600)

    @series begin
        subplot := 1
        yguide := "input"
        seriescolor --> SECONDARY_COLOR
        ns
    end
    @series begin
        subplot_offset := 1
        draws := draws
        ns.x, nm
    end
end


# -----------------------------------
# enable plotting of auto covariance
# -----------------------------------

"""
    plotautocovfit(ns::NoiseSample, nm::MvGaussianNoiseModel; kw...)

Plot results of autocovariance fit.
"""
@userplot PlotAutoCovFit
@recipe function plot_repice(pac::PlotAutoCovFit)

    ns, nm = pac.args

    xguide := "lag"
    yguide := "auto-covariance"

    lags, acov = estimate_autocov(ns)

    @series begin
        label := "estimate"
        seriescolor --> PRIMARY_COLOR
        lags, acov
    end

    @series begin
        label := @sprintf "fit (α = %.3e, λ = %.3e)" nm.α nm.λ
        seriescolor --> SECONDARY_COLOR
        lags, gauss_kernel(lags, [nm.α, nm.λ])
    end
end

MonteCarloMeasurements.mcplot(uc::UncertainCurve; draws=10, alpha=0.5, kw...) = MonteCarloMeasurements.mcplot(uc.x, uc.y, draws; alpha=0.5, kw...)


# --------------------------------------------
# enable plotting of UncertainBound histograms
# --------------------------------------------

@recipe plot_repice(::Type{T}, ub::T) where {T<:UncertainBound} = [ub.left.particles, ub.right.particles]

# --------------------------------------------
# enable plotting of Fits
# --------------------------------------------

@recipe function plot_recipe(f::PseudoVoigtFit{T}) where {T<:Number}
    w = WIDTH_FIT_PLOT
    #! format: off
    left        = f.center - f.width * w
    right       = f.center + f.width * w
    span        = abs(right - left)
    xs          = collect(left:(span/100):right)
    ys          = pvoigt_profile(xs, f)
    baseline    = xs .* f.slope .+ f.offset
    peak_height = pvoigt_peak(f.area, f.width, f.mixing)
    offset      = mean(baseline) + peak_height * OFFSET_FIT_PLOT
    #! format: on

    # fit
    @series begin
        color := :black
        label := nothing
        alpha := 0.5
        fill := (baseline, SECONDARY_COLOR)
        fillalpha := 0.2
        linewidth := 0.0
        xs, ys
    end
    # baseline corrected peak
    @series begin
        color := SECONDARY_COLOR
        label := nothing
        xs, ys .- baseline .+ offset
    end
    # line marking peak width
    @series begin
        color := SECONDARY_COLOR
        alpha := 0.25
        label := nothing
        [f.center - 0.5f.width, f.center + 0.5f.width], [1, 1] .* (offset + 0.5peak_height)
    end
    # line marking peak center and height
    @series begin
        color := SECONDARY_COLOR
        alpha := 0.25
        label := nothing
        [f.center, f.center], [offset, offset + peak_height]
    end
    # line marking baseline
    @series begin
        color := SECONDARY_COLOR
        alpha := 0.25
        label := nothing
        [f.center - w * f.width, f.center + w * f.width], [offset, offset]
    end
end

@recipe function plot_recipe(crv::Curve{T}, fs::Vector{PseudoVoigtFit{T}}) where {T<:Number}
    for f in fs
        @series begin
            f
        end
    end
    @series begin
        color --> PRIMARY_COLOR
        crv
    end

end

# plot draws of curves alongside with draws of fits
@recipe function plot_recipe(
    uc::UncertainCurve{T,N},
    fits::Vector{UncertainPseudoVoigtFit{T,N}}
    ;
    draws=3,
) where {T,N}

    # Be flexible here, if user gives a range, we draw these specific draws
    range = draws isa UnitRange ? draws : 1:draws

    legend := :none
    layout := (length(draws), 1)
    link := :both
    size --> (500, 600)

    for (k, i) ∈ enumerate(range)
        for (j, f) in enumerate(fits)
            @series begin
                fillcolor := j % 2 == 1 ? :red : :orange
                subplot := k
                bound := f
                get_draw(i, f)
            end
        end
        @series begin
            subplot := k
            seriescolor := PRIMARY_COLOR
            yguide := "sample $(i)"
            get_draw(i, uc)
        end
    end
end

# plot histograms of fit-parameters
@recipe function plot_recipe(
    fit::UncertainPseudoVoigtFit{T,N}
) where {T,N}

    alpha --> 0.5
    fill --> true
    layout --> @layout [
        a _ _ _ _ _;
        a a _ _ _ _;
        a a a _ _ _;
        a a a a _ _;
        a a a a a _;
        a a a a a a;
    ]
    legend --> false
    size --> (1000, 1000)
    xrotation := 60

    fields = (:area, :center, :width, :mixing, :offset, :slope)
    for (j, fj) in enumerate(fields)
        for (i, fi) in enumerate(fields)
            if j < i
                continue # we only want to plot the lower diagonal
            end
            if fi == fj
                @series begin
                    seriestype := :stephist
                    ylabel --> "counts"
                    if i == 6
                        xlabel := String(fi)
                    end
                    getfield(fit, fi)
                end
            else
                @series begin
                    seriestype := :scatter
                    if j == 6
                    xlabel := String(fi)
                    end
                    if i == 1
                        ylabel := String(fj)
                    end
                    getfield(fit, fi).particles, getfield(fit, fj).particles
                end
            end
        end
    end
end
