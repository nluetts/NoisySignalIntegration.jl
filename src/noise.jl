"""
This file defines the type `Noise` and basic functions
to analyse a noise sample.
"""

## helper functions

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

"""
    function get_cov(x::Array{Float64, 1}, α::Float64, λ::Float64)

Return covariance matrix `Σ` for the exponentiated quadratic kernel
with amplitude `α` and lag `λ`.
"""
function get_cov(x::Array{Float64, 1}, α::Float64, λ::Float64)
    N = length(x)
    Σ = Array{Float64}(undef, N, N)
    for i=1:N, j=1:N
        Σ[i,j] = (α^2)*exp(-0.5*((x[i] - x[j])/λ)^2) # exponentiated quadratic kernel
        if i == j
            Σ[i,j] *= 1.0000001 # for numerical stability (?)
        end
    end
    Σ
end

## datatypes to hold noise sample and parameters

"""
#TODO
"""
struct NoiseSample{T} <: AbstractSpectrum{T}
    
    x::Array{T, 1}
    y::Array{T, 1}
    lags::Array{T, 1}
    autocov::Array{T, 1}

    function NoiseSample(s::Spectrum{T}; detrend_order::Integer=0) where {T<:AbstractFloat}
        # check that spectral grid is equally spaced
        if !allapproxequal(diff(s.x))
            throw(ArgumentError("Noise analysis only works on an evenly spaced spectral grid."))
        end
        δx = abs(s.x[1]-s.x[2]) # x-values must be evenly spaced
        noise_sample_y = detrend_order > 0 ? detrend(s.x, s.y, detrend_order) : s.y
        lag_units = begin
            N = convert(Integer, round(length(s)/2))
            collect(1:N)
        end
        lags = δx .* lag_units
        acov = autocov(noise_sample_y, lag_units)
        return new{T}(s.x, noise_sample_y, lags, acov)
    end
end

function plot_autocov!(p::Plots.Plot, ns::NoiseSample, args...; kw...) 
    p = plot!(p, ns.lags, ns.autocov, args...; label="autocovariance", kw...)
    xlabel!(p, "lag")
    ylabel!(p, "autocovariance")
    return p
end
plot_autocov(ns::NoiseSample, args...; kw...) = plot_autocov!(plot(), ns, args...; kw...)

# fit function for noise autocorrelation
function gauss_kernel(Δx, β)
    α, λ = β
    return @. α^2 * exp(-0.5 * (Δx/λ)^2)
end

"""
    NoiseFit{T<:AbstractFloat}

Results of fitting autocovariance of a noise sample with a Gaussian kernel.

# Fields
- `sample::NoiseSample`
- `α::AbstractFloat`: amplitude
- `λ::AbstractFloat`: lag

# Constructor
    NoiseFit(ns::NoiseSample{T}; α_guess::T=1.0, λ_guess::T=1.0) where {T<:AbstractFloat}

Estimates parameters α and λ of the exponentiated quadratic covariance function
k(Δx) = α² exp(-0.5 (Δx/λ)²) by fitting the estimated autocovaraiance of the
noise sample `ns`.

Note: Change initial guess if fit does not converge to sensible result.
"""
struct NoiseFit{T<:AbstractFloat}
    sample::NoiseSample{T}
    α::T
    λ::T
    function NoiseFit(ns::NoiseSample{T}; α_guess=1.0, λ_guess=1.0) where {T}
        fit = curve_fit(gauss_kernel, ns.lags, ns.autocov, [α_guess, λ_guess])
        return new{T}(ns, fit.param...)
    end
end

function Base.show(io::IO, np::NoiseFit{T}) where {T}
    println("NoiseFit{$T} for sample with $(length(np.sample)) datapoints:")
    println("α = $(np.α), λ = $(np.λ)")
end

function Plots.plot!(p::Plots.Plot, np::NoiseFit, args...; kw...)
    fitlabel = @sprintf "fit (α = %.3e, λ = %.3e)" np.α np.λ
    plot_autocov!(p, np.sample, args...; kw...)
    fity = begin
        β = [np.α, np.λ]
        gauss_kernel(np.sample.lags, β)
    end
    plot!(p, np.sample.lags, fity, args...; kw..., label=fitlabel)
    return p
end
Plots.plot(np::NoiseFit) = plot!(plot(), np)

## functions to sample noise and plot result

function sample(
    np::NoiseFit,
    s::AbstractSpectrum;
    n::Int64=1,
    seed::Union{Nothing, Int64}=nothing
)
    !isnothing(seed) && seed!(seed)
    samples = rand(MvNormal(get_cov(s.x, np.α, np.λ)), n)
    n == 1 && return samples[:,1]
    return samples
end

function sample(np::NoiseFit; n::Int64=1, seed::Union{Nothing, Int64}=nothing)
    sample(np, np.sample; n=n, seed=seed)
end

function plot_samples(
    ns::NoiseFit,
    n_samples::Integer,
    args...;
    seed::Union{Nothing, Int64}=nothing,
    kw...
)
    n_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    p = plot(ns.sample, args..., label="experimental", kw...)
    y = ns.sample.y
    x = ns.sample.x
    span = maximum(y) - minimum(y)
    ss = sample(ns; n=n_samples, seed=seed)
    for (i, s) in enumerate(eachcol(ss))
        p = plot!(p, x, s .+ i*2*span, label="sample $(i)")
    end
    return p
end