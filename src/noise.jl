"""
    detrend(x, y, poly_order)

Subtract polynomial from `y` data.
"""
detrend(x, y, poly_order) = y - fit(x, y, poly_order).(x)
detrend(c::Curve, poly_order) = Curve(c.x, detrend(c.x, c.y, poly_order))


"""
    NoiseSample <: AbstractCurve

Curve holding a noise sample. At minimum, a constant offset
is removed from the noise sample upon construction.
If provided, a ploynomial of order `n` is subtracted.

# Fields
- `x::Vector{Float64}` : x-grid
- `y::Vector{Float64}` : y-value

# Constructors

    NoiseSample(x::Vector{T}, y::Vector{T})
    NoiseSample(x::Vector{T}, y::Vector{T}, n::Int)
    NoiseSample(s::Curve)
    NoiseSample(s::Curve, n::Int)

# Notes

Noise samples should not have a slow trend on top of the noise,
otherwise noise parameters cannot be reliably infered.
Construct `NoiseSample` objects with a high enough order
`n` so as to remove any slow variations.
"""
struct NoiseSample{T} <: AbstractCurve
    x::Vector{T}
    y::Vector{T}
    function NoiseSample{T}(x, y, n::Int=0) where T
        verify_same_length(x, y)
        return new{T}(x, detrend(x, y, n))
    end
end

NoiseSample(x::Vector{T}, y::Vector{T}, n::Int=0) where T = NoiseSample{T}(x, y, n)
NoiseSample(y::Vector{T}, n::Int=0) where T = NoiseSample(collect(T, 1:length(y)), y, n)
NoiseSample(s::Curve, n::Int=0) = NoiseSample(s.x, s.y, n)

# -------------------------------------
# Noise models
# -------------------------------------

"""
    AbstractNoiseModel

Supertype of noise models.
"""
abstract type AbstractNoiseModel end

"""
    GaussianNoiseModel{T <: Real} <: AbstractNoiseModel

Model to describe noise following a univariate Gaussian distribution
(uncorrelated noise).

# Fields

- `σ::T` : standard deviation
"""
struct GaussianNoiseModel{T <: Real} <: AbstractNoiseModel
    σ::T
end

Base.eltype(::GaussianNoiseModel{T}) where T = T

function Base.show(io::IO, nm::GaussianNoiseModel)
    println(io, "GaussianNoiseModel(σ = $(nm.σ))")
end

"""
    MvGaussianNoiseModel{T <: Real} <: AbstractNoiseModel

Model to describe noise following a multivariate Gaussian distribution
(correlated noise).

# Fields
- `δx::T`: noise grid spacing
- `α::T` : autocovariance amplitude
- `λ::T` : autocovaraiance lag
"""
struct MvGaussianNoiseModel{T <: Real} <: AbstractNoiseModel
    δx::T
    α::T
    λ::T
end
MvGaussianNoiseModel(δx, α, λ) = MvGaussianNoiseModel(promote(δx, α, λ)...)

Base.eltype(::MvGaussianNoiseModel{T}) where T = T

function Base.show(io::IO, nm::MvGaussianNoiseModel)
    println(io, "MvGaussianNoiseModel(α = $(nm.α), λ = $(nm.λ))")
end

"""
    gauss_kernel(Δx, β)

Fit function for noise autocorrelation.

    β = [α, λ]
    f(δx) = α^2 * exp(-0.5 * (Δx/λ)^2)
"""
function gauss_kernel(Δx, β)
    α, λ = β
    return @. α^2 * exp(-0.5 * (Δx/λ)^2)
end

"""
    get_cov(δx::T, n::Int, α::T, λ::T) where {T <: Real}

Return covariance matrix `Σ` with size `n` × `n` for the exponentiated quadratic kernel
with amplitude `α` and lag `λ` calculated on a grid of spacing `δx`.
"""
function get_cov(δx::T, n::Int, α::T, λ::T) where {T <: Real}
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
get_cov(nm::MvGaussianNoiseModel{T}, n::Int) where {T <: Real} = get_cov(nm.δx, n, nm.α, nm.λ)


function estimate_autocov(n::NoiseSample)
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

function _fit_noise(lags::Vector{T}, autocov::Vector{T}; α_guess=1.0, λ_guess=1.0) where {T <: Real}
    fit = curve_fit(gauss_kernel, lags, autocov, [α_guess, λ_guess])
    return MvGaussianNoiseModel(lags[1], fit.param...)
end

"""
    fit_noise(ns::NoiseSample; α_guess=1.0, λ_guess=1.0)

Build a ``MvGaussianNoiseModel` for the `NoiseSample` `ns`.

Estimates parameters `α` and `λ` of the exponentiated quadratic covariance function
k(Δx) = α² exp(-0.5 (Δx/λ)²) fitted to autocovariance (k) vs. lag (Δx) of the noise sample.

Change initial guess if fit does not converge to sensible result.
"""
function fit_noise(n::NoiseSample; α_guess=1.0, λ_guess=1.0)
    lags, acov = estimate_autocov(n)
    return _fit_noise(lags, acov; α_guess=α_guess, λ_guess=λ_guess)
end

"""
Generate correlated noise from a noise sample or noise model.
"""
function generate_noise(nm::MvGaussianNoiseModel{T}, len::Int, samples::Int) where {T <: Real}
    Σ = get_cov(nm, len)
    return Particles(samples, MvNormal(Σ))
end
function generate_noise(ns::NoiseSample, samples::Int; kw...)
    nm = fit_noise(ns; kw...)
    δx = ns.x[2] - ns.x[1]
    return generate_noise(nm, length(ns), samples)
end

"""
Generate uncorrelated noise from a GaussianNoiseModel.
"""
function generate_noise(nm::GaussianNoiseModel{T}, len::Int, samples::Int) where {T <: Real}
    return Particles(samples, MvNormal(zeros(T, len), nm.σ))
end


"""
Add noise from noise model to curve, retrieve UncertainCurve
"""
function add_noise(c::Curve, nm::GaussianNoiseModel, samples::Int=10_000)
    return UncertainCurve(c.x, Particles(samples, MvNormal(c.y, nm.σ)))
end
function add_noise(c::Curve, nm::MvGaussianNoiseModel, samples::Int=10_000)
    δx = c.x[2] - c.x[1]
    Σ = get_cov(nm, length(c))
    return UncertainCurve(c.x, Particles(samples, MvNormal(c.y, Σ)))
end


