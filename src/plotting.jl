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

@recipe function plot_recipe(nm::AbstractNoiseModel; grid_points::Integer=1000, noise_samples::Integer=3)
    noise_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    S = sample(nm, grid_points, noise_samples)
    span = (maximum(S) - minimum(S)) * 1.1
    for (i, s) in enumerate(eachcol(S))
        @series begin
            label := "sample $(i)"
            s .+ (i*span - 1)
        end
    end
end

@recipe function plot_recipe(ns::Noise, nm::AbstractNoiseModel; noise_samples::Integer=3)
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