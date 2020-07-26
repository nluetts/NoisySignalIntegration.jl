# -------------------------------------
# Curves
# -------------------------------------

"""
    AbstractCurve

Abstract supertype for all curves (xy-data).
Subtypes need to implement fields `x` and `y`.
"""
abstract type AbstractCurve end

Base.length(s::AbstractCurve) = length(s.x)
Base.lastindex(s::AbstractCurve) = length(s)
Base.getindex(s::AbstractCurve, i::Integer) = (s.x[i], s.y[i])
Base.getindex(s::AbstractCurve, r::UnitRange) =Curve(s.x[r], s.y[r])
Base.iterate(s::AbstractCurve, i::Int64) = i > length(s) ? nothing : (s[i], i + 1)
Base.iterate(s::AbstractCurve) = iterate(s, 1)
Base.IteratorSize(itertype::Type{AbstractCurve}) = Base.HasLength()
Base.minimum(s::AbstractCurve) = s[argmin(s.y)]
Base.maximum(s::AbstractCurve) = s[argmax(s.y)]

Base.:(==)(c0::AbstractCurve, c1::AbstractCurve) = c0.x == c1.x && c0.y == c1.y
Base.:(+)(c::AbstractCurve, y) = typeof(c)(c.x, c.y .+ Float64.(y))
Base.vcat(c0::AbstractCurve, c1::AbstractCurve) = typeof(c0)(vcat(c0.x, c1.x), vcat(c0.y, c1.y))

function Base.show(io::IO, c::AbstractCurve)
    println("$(typeof(c)), $(length(c)) datapoints")
    if length(c) > 10
        ps = [c[1:5]; c[end-4:end]]
        indicator = "    ⋮\n"
    else
        ps = c
        indicator = ""
    end
    for (i, p) in enumerate(ps)
        @printf "    %-5e %-5e\n" p...
        i == 5 ? print(indicator) : nothing
    end
end

"""
    Curve <: AbstractCurve

Datatype holding x-y data.
x and y vectors have to have the same length.

# Fields
- `x::Vector{Float64}` : spectral grid
- `y::Vector{Float64}` : spectral intensity

# Constructors

    Curve(x::Vector{Float64}, y::Vector{Float64})
"""
struct Curve <: AbstractCurve
    x::Vector{Float64}
    y::Vector{Float64}
    function Curve(x::Vector{Float64}, y::Vector{Float64})
        verify_same_length(x, y)
        return new(x, y)
    end
end

"""
    Noise <: AbstractCurve

Curve holding a noise sample. At minimum, a constant offset
is removed from the noise sample upon construction.
If provided, a ploynomial of order `n` is subtracted.

# Fields
- `x::Vector{Float64}` : spectral grid
- `y::Vector{Float64}` : spectral intensity

# Constructors
    Noise(x::Vector{Float64}, y::Vector{Float64})
    Noise(x::Vector{Float64}, y::Vector{Float64}, n::Integer)
    Noise(s::Curve)
    Noise(s::Curve, n::Integer)

"""
struct Noise <: AbstractCurve
    x::Vector{Float64}
    y::Vector{Float64}
    function Noise(x, y)
        verify_same_length(x, y)
        return new(x, detrend(Float64.(x), Float64.(y), 0))
    end
end

Noise(x, y, n::Integer) = Noise(x, detrend(x, y, n))
Noise(s::Curve, n::Integer=0) = Noise(s.x, s.y, n)

# -------------------------------------
# Noise models
# -------------------------------------

"""
    AbstractNoiseModel

Supertype of noise models.
Subtypes must support the following methods:

- `sample(nm::NoiseModel, m::Integer [, n::Integer])` :: Array{Float64, 2}
  with:
  `m` = no. of data points per noise sample
  `n` = no. of noise samples
"""
abstract type AbstractNoiseModel end

"""
    GaussianNoiseModel

Model to describe noise following a univariate Gaussian distribution
(uncorrelated noise).

# Fields
- `σ::T` : standard deviation
"""
struct GaussianNoiseModel <: AbstractNoiseModel
    σ::Float64
    GaussianNoiseModel(σ) = new(Float64(σ))
end
function Base.show(io::IO, nm::GaussianNoiseModel)
    println("GaussianNoiseModel(σ = $(nm.σ))")
end

"""
    MvGaussianNoiseModel

Model to describe noise following a multivariate Gaussian distribution
(correlated noise).

# Fields
- `δx::Float64`: noise grid spacing
- `α::Float64` : autocovariance amplitude
- `λ::Float64` : autocovaraiance lag
"""
struct MvGaussianNoiseModel <: AbstractNoiseModel
    δx::Float64
    α::Float64
    λ::Float64
end

function Base.show(io::IO, nm::MvGaussianNoiseModel)
    println("MvGaussianNoiseModel(α = $(nm.α), λ = $(nm.λ))")
end

# -------------------------------------
# Bounds
# -------------------------------------

"""
    ScaledShiftedBeta

Scaled and shifted `Beta(α, β)` distribution.
Samples fall in the interval [`a`, `b`].
"""
ScaledShiftedBeta = LocationScale{Float64, Beta{Float64}}

"""
    scale_shift_beta(α::Float64, β::Float64, a::Float64, b::Float64)

Create a scaled and shifted `Beta(α, β)` distribution.
Samples fall in the interval [`a`, `b`].
"""
function scale_shift_beta(α, β, a, b)
    return LocationScale(a, b - a, Beta(α, β))
end

"""
    AbstractUncertainBound

Supertype of uncertain bounds.
Subtypes must implement methods to sample integration bounds.
"""
abstract type AbstractUncertainBound end

struct LeftRightBound{S1<:ContinuousUnivariateDistribution,
                      S2<:ContinuousUnivariateDistribution} <: AbstractUncertainBound
    left::S1
    right::S2
end

# sample cache that is retrieved by `WidthBoundClone`s
const WIDTH_SAMPLES = Dict{Int, Array{Float64}}()

struct WidthBound{T<:ContinuousUnivariateDistribution} <: AbstractUncertainBound
    loc::Float64
    width::T
    id::Int
    function WidthBound(loc::Float64, width::T) where {T<:ContinuousUnivariateDistribution}
        id = isempty(WIDTH_SAMPLES) ? 1 : maximum(keys(WIDTH_SAMPLES)) + 1
        WIDTH_SAMPLES[id] = []
        return new{T}(loc, width, id)
    end
end

struct WidthBoundClone <: AbstractUncertainBound
    loc::Float64
    reference::WidthBound
end

WidthBoundUnion = Union{WidthBound, WidthBoundClone}