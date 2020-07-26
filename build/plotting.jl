# --------------------------------------------
# enable plotting of custom types via recipes
# --------------------------------------------


# --------------------------------------------
# enable plotting of curves and noise samples
# --------------------------------------------
@recipe function plot_recipe(crv::AbstractCurve)
    xguide --> "x"
    yguide --> "y"
    return crv.x, crv.y
end


# --------------------------------------
# enable plotting of noise sample draws
# --------------------------------------

@recipe function plot_recipe(nm::AbstractNoiseModel; grid_points=1000, noise_samples=3)
    noise_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    S = sample(nm, grid_points, noise_samples)
    delete!(plotattributes, :grid_points)
    delete!(plotattributes, :noise_samples)
    span = (maximum(S) - minimum(S)) * 1.1
    for (i, s) in enumerate(eachcol(S))
        @series begin
            label := "sample $(i)"
            s .+ (i*span - 1)
        end
    end
end

@recipe function plot_recipe(ns::Noise, nm::AbstractNoiseModel; noise_samples=3)
    noise_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    S = sample(nm, length(ns), noise_samples)
    span = (maximum(S) - minimum(S)) * 1.1
    for (i, s) in enumerate(eachcol(S))
        @series begin
            label := "sample $(i)"
            ns.x, s .+ i*span
        end
    end
    seriescolor --> :black
    label := "experimental"
    ns.x, ns.y
end

# -----------------------------------
# enable plotting of auto covariance
# -----------------------------------
struct AutoCovPlot end # just for dispatch

@recipe function plot_repice(ns::Noise, nm::MvGaussianNoiseModel, ::Type{AutoCovPlot})
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

plot_autocov(ns::Noise, nm::MvGaussianNoiseModel; kw...) = plot(ns, nm, AutoCovPlot; kw...)

# --------------------------------------
# enable plotting of integration bounds
# --------------------------------------

@recipe function plot_recipe(crv::Curve, left::Float64, right::Float64)
    fillrange := 0
    fillalpha --> 0.5
    fillcolor --> :orange
    linewidth --> 0
    label     --> nothing

    bounds_parameters = get_integration_bounds(crv.x, crv.y, left, right)
    l, r = bounds_parameters[2:3]
    x_l, x_r, y_l, y_r = bounds_parameters[12:15]
    x = [x_l; crv.x[l+1:r-1]; x_r; x_l]
    y = [y_l; crv.y[l+1:r-1]; y_r; y_l]
    return x, y
end

@recipe function plot_recipe(crv::Curve, bnds::Vector{T}) where {T <: AbstractUncertainBound}
    for b in bnds
        @series begin
            l, r = typeof(b) <: WidthBoundUnion ? sample(b, crv) : sample(b)
            crv, l, r
        end
    end
    return crv
end

@recipe function plot_recipe(crv::Curve, bnds::Vector{T}, nm::AbstractNoiseModel; samples=3) where {T <: AbstractUncertainBound}
    
    legend --> :outertopright

    _, y_min = minimum(crv)
    _, y_max = maximum(crv)

    span = (y_max - y_min) * 1.1

    # plot experimental spectrum with mean bounds
    @series begin
        label := "experiment, mean bounds"
        crv
    end
    for b in bnds
        @series begin
            if typeof(b) <: WidthBound
                l, r, = left_right_from_peak(crv.x, crv.y, b.loc, mean(b.width))
            elseif typeof(b) == WidthBoundClone
                :fillcolor --> :red
                l, r, = left_right_from_peak(crv.x, crv.y, b.loc, mean(b.reference.width))
            else
                l, r = map(mean, [b.left, b.right])
            end
            crv, l, r
        end
    end

    # plot samples
    for i in 1:samples
        label --> "sample $i"
        crv_ = crv + (sample(nm, length(crv)) .+ span * i)
        @series begin
            crv_
        end
    
        for b in bnds
            @series begin
                label := nothing
                if typeof(b) == WidthBoundClone
                    :fillcolor --> :red
                end
                l, r = typeof(b) <: WidthBoundUnion ? sample(b, crv_) : sample(b)
                crv_, l, r
            end
        end
    end
end