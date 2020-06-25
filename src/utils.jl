"""
    crop(s::Curve{T}, left::T, right::T) where {T}

Crop a slice from a spectrum and return a new spectrum.
"""
function crop(s::Curve{T}, left::T, right::T) where {T}
    i = searchsortedlast(s.x, left)
    j = searchsortedfirst(s.x, right)
    return Curve(s.x[i:j], s.y[i:j])
end

Plots.plot!(p::Plots.Plot, s::AbstractCurve, args...; kw...) = plot!(p, s.x, s.y, args...; kw...)
Plots.plot(s::AbstractCurve, args...; kw...) = plot!(plot(), s.x, s.y, args...; kw...)

#######################################################
### plot autocovariance ###############################
#######################################################

function plot_autocov!(
    p::Plots.Plot,
    lags::Vector{T},
    autocov::Vector{T},
    args...;
    kw...
) where {T<:AbstractFloat}
    p = plot!(p, lags, autocov, args...; label="autocovariance", kw...)
    xlabel!(p, "lag")
    ylabel!(p, "autocovariance")
    return p
end
function plot_autocov(
    lags::Vector{T},
    autocov::Vector{T},
    args...;
    kw...
) where {T<:AbstractFloat}
    return plot_autocov!(plot(), lags, autocov, args...; kw...)
end

function plot_autocov!(p::Plots.Plot, n::Noise, args...; kw...)
    lags, acov = estimate_autocov(n)
    return plot_autocov!(p, lags, acov, args...; kw...)
end
plot_autocov(n::Noise, args...; kw...) = plot_autocov!(plot(), n, args...; kw...)

function plot_autocov!(p::Plots.Plot, n::Noise, nm::MvGaussianNoiseModel, args...; kw...)
    lags, acov = estimate_autocov(n)
    p = plot_autocov!(p, lags, acov, args...; kw...)
    fity = begin
        β = [nm.α, nm.λ]
        gauss_kernel(lags, β)
    end
    fitlabel = @sprintf "fit (α = %.3e, λ = %.3e)" nm.α nm.λ
    return plot!(p, lags, fity, args...; kw..., label=fitlabel)
end
function plot_autocov(n::Noise, nm::MvGaussianNoiseModel, args...; kw...)
    return plot_autocov!(plot(), n, nm, args...; kw...)
end

#######################################################
### plot noise samples ################################
#######################################################


function Plots.plot!(p::Plots.Plot, nm::AbstractNoiseModel; grid_points::Integer=1000, noise_samples::Integer=3)
    noise_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    S = sample(nm, grid_points, noise_samples)
    span = (maximum(S) - minimum(S)) * 1.1
    for (i, s) in enumerate(eachcol(S))
        p = plot!(p, s .+ i*span, label="sample $(i)")
    end
    return p
end
function Plots.plot(nm::AbstractNoiseModel; grid_points::Integer=100, noise_samples::Integer=1)
    return plot!(plot(), nm; grid_points=grid_points, noise_samples=noise_samples)
end

function Plots.plot!(
    p::Plots.Plot,
    n::Noise,
    nm::AbstractNoiseModel,
    args...;
    noise_samples::Integer=3,
    kw...
)
    noise_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    S = sample(nm, length(n.x), noise_samples)
    span = (maximum(S) - minimum(S)) * 1.1
    for (i, s) in enumerate(eachcol(S))
        p = plot!(p, n.x, s .+ i*span, label="sample $(i)")
    end
    offset = span * (size(S, 2) + 1)
    return plot!(p, n.x, n.y .+ offset, args..., label="experimental", kw...)
end
function Plots.plot(
    n::Noise,
    nm::AbstractNoiseModel,
    args...;
    noise_samples::Integer=3,
    kw...
)
    return plot!(plot(), n, nm, args...; noise_samples=noise_samples, kw...)
end

#######################################################
### plot bounds #######################################
#######################################################

function _plot_bound!(p::Plots.Plot, crv::Curve{T}, left::T, right::T) where {T}
    bounds_parameters = get_integration_bounds(crv.x, crv.y, left, right)
    l, r = bounds_parameters[2:3]
    x_l, x_r, y_l, y_r = bounds_parameters[12:15]
    x_ = [x_l; crv.x[l+1:r-1]; x_r; x_l]
    y_ = [y_l; crv.y[l+1:r-1]; y_r; y_l]
    return plot!(p, x_, y_; fill=(0, 0.5, :orange), linewidth=0)
end

