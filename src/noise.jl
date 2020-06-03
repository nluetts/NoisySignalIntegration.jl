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
struct NoiseSample{T<:AbstractFloat}
    
    spectrum::Spectrum{T}
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
        return new{T}(Spectrum(s.x, noise_sample_y), lags, acov)
    end
end

Base.length(n::NoiseSample) = length(n.spectrum)

function Base.show(io::IO, n::NoiseSample{T}) where {T}
    println("NoiseSample holding")
    show(io, n.spectrum) 
end

Plots.plot!(p::Plots.Plot, ns::NoiseSample, args...; kw...) = plot!(p, ns.spectrum, args...; kw...)
Plots.plot(ns::NoiseSample, args...; kw...) = plot!(plot(), ns, args...; kw...)


"""
    NoiseParameters(sample::NoiseSample{T}, α::T, λ::T) where {T<:AbstractFloat}

Datatype holding a noise sample and the noise amplitude `α`
and autocovariance lag `λ`.
"""
struct NoiseParameters{T<:AbstractFloat}
    sample::NoiseSample{T}
    α::T
    λ::T
end

function Base.show(io::IO, np::NoiseParameters{T}) where {T}
    println("NoiseParameters{$T} for sample with $(length(np.sample.spectrum)) datapoints:")
    println("α = $(np.α), λ = $(np.λ)")
end

# fit function for noise autocorrelation
function gauss_kernel(Δx, β)
    α, λ = β
    return @. α^2 * exp(-0.5 * (Δx/λ)^2)
end

"""
    fit_noise(ns::NoiseSample; α_guess=1.0, λ_guess=1.0) :: NoiseParameters

Estimate parameters α and λ of the exponentiated quadratic covariance function
k(Δx) = α² exp(-0.5 (Δx/λ)²) by fitting the estimated autocovaraiance of the
noise sample `ns`.

Note: Change initial guess if fit does not converge to sensible result.
"""
function fit_noise(ns::NoiseSample; α_guess=1.0, λ_guess=1.0) :: NoiseParameters
    fit = curve_fit(gauss_kernel, ns.lags, ns.autocov, [α_guess, λ_guess])
    return NoiseParameters(ns, fit.param...)
end

function Plots.plot!(p::Plots.Plot, np::NoiseParameters, args...; kw...)
    kw_exp = merge(Dict(:label => "autocovariance"), kw)
    fitlabel = @sprintf "fit (α = %.3e, λ = %.3e)" np.α np.λ
    kw_fit = merge(kw, Dict(:label => fitlabel))
    plot!(p, np.sample.lags, np.sample.autocov, args...; kw_exp...)
    β = [np.α, np.λ]
    plot!(p, np.sample.lags, gauss_kernel(np.sample.lags, β), args...; kw_fit...)
    xlabel!(p, "lag")
    ylabel!(p, "autocovariance")
    return p
end
Plots.plot(np::NoiseParameters) = plot!(plot(), np)

## functions to sample noise and plot result

function sample(
    np::NoiseParameters,
    s::Spectrum;
    n::Int64=1,
    seed::Union{Nothing, Int64}=nothing
)
    !isnothing(seed) && seed!(seed)
    samples = rand(MvNormal(get_cov(s.x, np.α, np.λ)), n)
    n == 1 && return samples[:,1]
    return samples
end

function sample(np::NoiseParameters; n::Int64=1, seed::Union{Nothing, Int64}=nothing)
    sample(np, np.sample.spectrum; n=n, seed=seed)
end

function plot_samples(
    ns::NoiseParameters,
    n_samples::Integer,
    args...;
    seed::Union{Nothing, Int64}=nothing,
    kw...
)
    n_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    p = plot(ns.sample, args..., label="experimental", kw...)
    y = ns.sample.spectrum.y
    x = ns.sample.spectrum.x
    span = maximum(y) - minimum(y)
    ss = sample(ns; n=n_samples, seed=seed)
    for (i, s) in enumerate(eachcol(ss))
        p = plot!(p, x, s .+ i*2*span, label="sample $(i)")
    end
    return p
end