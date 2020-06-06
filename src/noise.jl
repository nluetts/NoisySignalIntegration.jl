"""
This file defines the type `Noise` and basic functions
to analyse a noise sample.
"""

#######################################################
### helper functions ##################################
#######################################################

"""
Check if all elements in `x` are approximately equal.
"""
function allapproxequal(x)
    length(x) < 2 && return true
    e1 = x[1]
    for i = 2:length(x)
        x[i] ≈ e1 || return false
    end
    return true
end

"""
    detrend(x, y, poly_order)

Subtract polynomial baseline from y data.
"""
detrend(x, y, poly_order) = y - fit(x, y, poly_order).(x)

# fit function for noise autocorrelation
function gauss_kernel(Δx, β)
    α, λ = β
    return @. α^2 * exp(-0.5 * (Δx/λ)^2)
end

"""
    get_cov(δx::T, n::Integer, α::T, λ::T) where {T<:AbstractFloat}

Return covariance matrix `Σ` with size `n` × `n` for the exponentiated quadratic kernel
with amplitude `α` and lag `λ` calculated on a grid of spacing `δx`.
"""
function get_cov(δx::T, n::Integer, α::T, λ::T) where {T<:AbstractFloat}
    β = [α, λ]
    Σ = Array{T}(undef, n, n)
    for i = (1:n), j = (1:n)
        Δx = δx * (i - j)
        Σ[i,j] = gauss_kernel(Δx, β)
        if i == j
            Σ[i,j] *= 1.0000001 # for numerical stability (?)
        end
    end
    return Σ
end

#######################################################
### definition and sampling of noise models ###########
#######################################################


"""
    AbstractNoiseModel{T<:AbstractFloat}

Supertype of noise models.
Subtypes must support the following methods:

- `sample(nm::NoiseModel{T}, m::Integer [, n::Integer])` :: Array{T, N}
  with:
  `m` = no. of data points per noise sample
  `n` = no. of noise samples
"""
abstract type AbstractNoiseModel{T<:AbstractFloat} end

"""
    GaussianNoiseModel{T<:AbstractFloat}

Model to describe noise following a univariate Gaussian distribution
(uncorrelated noise).

# Fields
- `σ::T` : standard deviation
"""
struct GaussianNoiseModel{T} <: AbstractNoiseModel{T}
    σ::T
end
function Base.show(io::IO, nm::GaussianNoiseModel{T}) where {T}
    print("GaussianNoiseModel{$T}(σ = $(nm.σ))")
end

function sample(nm::GaussianNoiseModel{T}, m::Integer=100, n::Integer=1) where {T}
    return nm.σ * randn(T, m, n)
end

"""
    MvGaussianNoiseModel{T<:AbstractFloat}

Model to describe noise following a multivariate Gaussian distribution
(correlated noise).

# Fields
- `δx::T`: noise grid spacing
- `α::T` : autocovariance amplitude
- `λ::T` : autocovaraiance lag
"""
struct MvGaussianNoiseModel{T} <: AbstractNoiseModel{T}
    δx::T
    α::T
    λ::T
end

function Base.show(io::IO, nm::MvGaussianNoiseModel{T}) where {T}
    print("MvGaussianNoiseModel{$T}(α = $(nm.α), λ = $(nm.λ))")
end

function sample(
    nm::MvGaussianNoiseModel{T},
    m::Integer=100,
    n::Integer=1
) where {T<:AbstractFloat}
    samples = rand(MvNormal(get_cov(nm.δx, m, nm.α, nm.λ)), n)
    n == 1 && return samples[:,1]
    return samples
end

#######################################################
### fitting and plotting of noise sample ##############
#######################################################

"""
    Noise{T} <: AbstractSpectrum{T}

Spectrum holding a noise sample.

# Fields
- `x::T` : spectral grid
- `y::T` : spectral intensity

# Constructor
    Noise(s::Spectrum{T}; detrend_order::Integer=0) where {T<:AbstractFloat}

Fits a polynomial of order `detrend_order` to remove slow variation of noise
and returns a `Noise` instance holding the detrended data.
"""
struct Noise{T} <: AbstractSpectrum{T}
    x::Vector{T}
    y::Vector{T}
    function Noise(s::Spectrum{T}; detrend_order::Integer=0) where {T<:AbstractFloat}
        return new{T}(s.x, detrend(s.x, s.y, detrend_order))
    end
end

function estimate_autocov(n::Noise{T}) where {T<:AbstractFloat}
    # check that spectral grid is equally spaced
    if !allapproxequal(diff(n.x))
        throw(ArgumentError("Noise analysis only works on an evenly spaced spectral grid."))
    end
    δx = abs(n.x[1]-n.x[2]) # x-values must be evenly spaced
    lag_units = begin
        N = convert(Integer, round(length(n)/2))
        collect(1:N)
    end
    lags = δx .* lag_units
    acov = autocov(n.y, lag_units)
    return lags, acov
end

# fitting
"""
    fit_noise(lags::T, autocov::T; α_guess=1.0, λ_guess=1.0) where {T<:AbstractFloat}
    fit_noise(n::Noise{T}; α_guess=1.0, λ_guess=1.0) where {T}

Estimates parameters α and λ of the exponentiated quadratic covariance function
k(Δx) = α² exp(-0.5 (Δx/λ)²) fitted to autocovariance (k) vs. lag (Δx).
If a `Noise` instance is passes as first argument, its autocovarince is
estimated and passed to the fitting function.

Note: Change initial guess if fit does not converge to sensible result.
"""
function fit_noise(
    lags::Vector{T},
    autocov::Vector{T};
    α_guess=1.0,
    λ_guess=1.0
) where {T<:AbstractFloat}
    fit = curve_fit(gauss_kernel, lags, autocov, [α_guess, λ_guess])
    return MvGaussianNoiseModel{T}(lags[1], fit.param...)
end

function fit_noise(n::Noise{T}; α_guess=1.0, λ_guess=1.0) where {T}
    lags, acov = estimate_autocov(n)
    fit_noise(lags, acov; α_guess=α_guess, λ_guess=λ_guess)
end

# plot autocovariance

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

# plot noise samples

function Plots.plot!(p::Plots.Plot, nm::AbstractNoiseModel; grid_points::Integer=100, noise_samples::Integer=3)
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

