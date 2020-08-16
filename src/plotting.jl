const DEFAULT_COLOR = palette(:default)[1]


function get_left_right_points(x::AbstractArray{T}, y::AbstractArray{T}, left::T, right::T; subtract_baseline=true) where {T<:AbstractFloat}

    # this function is required to plot integration areas

    left, right = left < right ? (left, right) : (right, left)

    N = length(x)
    N != length(y) && throw(ArgumentError("Data arrays `x` and `y` must have the same length."))
    N < 2 && throw(ArgumentError("`x`and `y` must each hold at least two elements."))
    
    l =  r = xl = xr = yl = yr = nothing
    j = 2

    while j <= N

        x₀ = x[j-1]
        x₁ = x[j]
        y₀ = y[j-1]
        y₁ = y[j]
        
        if x₁ <= left
            j += 1
            continue
        elseif yl === nothing
            # this will only run once, when we enter the integration window
            if x₀ < left
                xl = left
                yl = lininterp(left, x₀, x₁, y₀, y₁)
                l = j - 1
            elseif x₀ == left
                xl = left
                yl = y₀
                l = j
            else
                # this case means that `left` <= x[1]
                xl = x₀
                yl = y₀
                l = 0
            end
        end
        
        if x₁ >= right
            # we move out of the integration window
            yr = x₁ == right ? y₁ : lininterp(right, x₀, x₁, y₀, y₁)
            xr = right
            r = j + 1
            break
        end
        j += 1
    end

    if yr === nothing
        # this case means that right > x[end]
        xr = x[N]
        yr = y[N]
        r = N + 1
    end
    return l, r, xl, xr, yl, yr
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
    layout := @layout grid(draws + 1, 1)
    link := :both
    
    delete!(plotattributes, :draws)
  
    for i ∈ 0:draws
        @series begin
            subplot := i + 1
            if i == 0
                seriescolor := :red
                yguide := "input"
                c
            else
                yguide := "sample $(i)"
                get_draw(i, uc)
            end
        end
    end
end
 

@recipe function plot_recipe(crv::Curve{T}, left::T, right::T; subtract_baseline=true) where T
    fillrange := 0
    fillalpha --> 0.5
    fillcolor --> :orange
    linewidth --> 0
    label     --> nothing
    
    left = T(left)
    right = T(right)
    
    l, r, xl, xr, yl, yr = get_left_right_points(crv.x, crv.y, left, right)
    if subtract_baseline
        x = [xl; crv.x[l+1:r-1]; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; yl]
    else
        x = [xl; crv.x[l+1:r-1]; xr; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; 0;  0 ]
    end
    return x, y
end


@recipe function plot_recipe(crv::Curve{T}, bnds::Vector{UncertainBound{T, N}}, draw::Int; subtract_baseline=true) where {T, N}
    @series begin
        crv
    end
    for b ∈ bnds
        @series begin
            subtract_baseline := subtract_baseline
            left, right = get_draw(draw, b)
            crv, left, right
        end
    end
end


@recipe function plot_recipe(crv::Curve{T}, bnd::UncertainBound{T, N}, draw::Int; subtract_baseline=true) where {T, N}
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
    layout := @layout grid(draws + 1, 1)
    link := :both
    
    mean_uc = mean(uc)
    
    for i ∈ 0:draws
        for b in bnds
            @series begin
                subplot := i + 1
                subtract_baseline := subtract_baseline
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
                seriescolor := :red
                yguide := "mean"
                mean_uc
            else
                seriescolor := DEFAULT_COLOR
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
    draws=3,
    subtract_baseline=true
) where {T, N}
    draws := draws
    subtract_baseline := subtract_baseline
    uc, [bnd]
end



# --------------------------------------
# enable plotting of noise sample draws
# --------------------------------------

@recipe function plot_recipe(x::Vector{T}, nm::AbstractNoiseModel; draws=3, subplot_offset=0) where {T <: Real}
    draws < 0 && throw(ArgumentError("Number of samples must be > 0."))

    layout --> @layout grid(draws, 1)
    legend --> :none
    link --> :both

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
    
    layout --> @layout grid(draws+1, 1)
    legend --> :none
    link --> :both

    @series begin
            subplot := 1
            yguide := "input"
            seriescolor --> :red
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
struct _AutoCovPlot end # just for dispatch

@recipe function plot_repice(ns::NoiseSample, nm::MvGaussianNoiseModel, ::Type{_AutoCovPlot})
    xguide := "lag"
    yguide := "auto-covariance"

    lags, acov = estimate_autocov(ns)

    @series begin
        label := "estimate"
        seriescolor --> :blue
        lags, acov
    end

    @series begin
        label := @sprintf "fit (α = %.3e, λ = %.3e)" nm.α nm.λ
        seriescolor --> :red
        lags, gauss_kernel(lags, [nm.α, nm.λ])
    end
end

plot_autocov(ns::NoiseSample, nm::MvGaussianNoiseModel; kw...) = plot(ns, nm, _AutoCovPlot; kw...)


MonteCarloMeasurements.mcplot(uc::UncertainCurve; draws=10, alpha=0.5, kw...) = MonteCarloMeasurements.mcplot(uc.x, uc.y, draws; alpha=0.5, kw...)