const PRIMARY_COLOR = :blue
const SECONDARY_COLOR = :red


# this function is required to plot integration areas
function get_left_right_points(
    xs::AbstractArray{T},
    ys::AbstractArray{T},
    xₗ::T,
    xᵣ::T,
    b::UncertainBound
    ;
    baseline_handeling=nothing
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
    ) && throw(error("At least one integration bound is outside the support region ($(minimum(x)), $(maximum(x)))."))


    if baseline_handeling == "local"
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
 

@recipe function plot_recipe(crv::Curve{T}, left::T, right::T; subtract_baseline=false, local_baseline=false, bound=nothing) where T
    fillrange := 0
    fillalpha --> 0.5
    fillcolor --> :orange
    linewidth --> 0
    label     --> nothing
    
    (local_baseline && bound == nothing) && error("You have to provide a bound if local_baseline == true.") |> throw
    (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

    left = T(left)
    right = T(right)
    
    if subtract_baseline
        bh = "end-to-end"
    elseif local_baseline
        bh = "local"
    else
        bh = nothing
    end
    
    l, r, xl, xr, yl, yr = get_left_right_points(crv.x, crv.y, left, right, bound; baseline_handeling=bh)
    x = [xl; crv.x[l+1:r-1]; xr; xl]
    if !(local_baseline || subtract_baseline)
        y = [yl; crv.y[l+1:r-1]; zero(T); zero(T)]
    else
        y = [yl; crv.y[l+1:r-1]; yr; yl]
    end
    
    return x, y
end


@recipe function plot_recipe(crv::Curve{T}, bnds::Vector{UncertainBound{T, N}}, draw::Int) where {T, N}
    @series begin
        crv
    end
    for b ∈ bnds
        @series begin
            left, right = get_draw(draw, b)
            bound := b
            crv, left, right
        end
    end
end


@recipe function plot_recipe(crv::Curve{T}, bnd::UncertainBound{T, N}, draw::Int) where {T, N}
    @series begin
        crv, [bnd], draw
    end
end


# plot draws of curves alongside with draws of bounds
@recipe function plot_recipe(
    uc::UncertainCurve{T, N},
    bnds::Vector{UncertainBound{T, N}}
    ;
    draws=3,
    subtract_baseline=true
) where {T, N}

    legend := :none
    layout := (draws + 1, 1)
    link := :both
    size --> (500, 600)
    
    mean_uc = mean(uc)
    
    for i ∈ 0:draws
        for (j, b) in enumerate(bnds)
            @series begin
                fillcolor := j % 2 == 1 ? :red : :orange
                subplot := i + 1
                bound := b
                if i == 0
                    mean_uc, mean(b)...
                else
                    get_draw(i, uc), get_draw(i, b)...
                end
            end
        end
        @series begin
            subplot := i + 1
            if i == 0
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
    uc::UncertainCurve{T, N},
    bnd::UncertainBound{T, N}
    ;
    draws=3
) where {T, N}
    draws := draws
    uc, [bnd]
end



# --------------------------------------
# enable plotting of noise sample draws
# --------------------------------------

@recipe function plot_recipe(x::Vector{T}, nm::AbstractNoiseModel; draws=3, subplot_offset=0) where {T <: Real}
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
    
    layout --> (draws+1, 1)
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

@recipe plot_repice(::Type{T}, ub::T) where {T <: UncertainBound} = [ub.left.particles, ub.right.particles]