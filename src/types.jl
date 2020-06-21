#######################################################
### Curves ############################################
#######################################################

"""
    AbstractCurve{T<:AbstractFloat}

Abstract supertype for all curves (xy-data).
Subtypes need to implement fields `x` and `y`.
"""
abstract type AbstractCurve{T<:AbstractFloat} end

Base.length(s::AbstractCurve) = length(s.x)
Base.lastindex(s::AbstractCurve) = length(s)
Base.getindex(s::AbstractCurve, i::Integer) = (s.x[i], s.y[i])
Base.getindex(s::AbstractCurve, r::UnitRange) =Curve(s.x[r], s.y[r])
Base.iterate(s::AbstractCurve, i::Int64) = i > length(s) ? nothing : (s[i], i + 1)
Base.iterate(s::AbstractCurve) = iterate(s, 1)
Base.IteratorSize(itertype::Type{AbstractCurve}) = Base.HasLength()

Base.:(==)(c0::AbstractCurve, c1::AbstractCurve) = c0.x == c1.x && c0.y == c1.y
Base.:(+)(c::AbstractCurve{T}, y::T) where {T} = typeof(c)(c.x, c.y + y)
Base.:(+)(c::AbstractCurve{T}, y::Vector{T}) where {T} = typeof(c)(c.x, c.y .+ y)
Base.vcat(c0::AbstractCurve, c1::AbstractCurve) = typeof(c0)(vcat(c0.x, c1.x), vcat(c0.y, c1.y))

function Base.show(io::IO, c::AbstractCurve{T}) where {T}
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
    Curve(x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
"""
struct Curve{T} <: AbstractCurve{T}
    x::Vector{T}
    y::Vector{T}
    function Curve{T}(x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
        verify_same_length(x, y)
        return new{T}(x, y)
    end
end

Curve(x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat} = Curve{T}(x, y)

"""
    Noise{T} <: AbstractCurve{T}

Curve holding a noise sample. At minimum, a constant offset is removed
from the noise sample upon construction.

# Fields
- `x::T` : spectral grid
- `y::T` : spectral intensity
- `n::Integer` (optional) : order of polynomial to subtract to detrend
  (defaults to `0`, i.e. constant offset removal)

# Constructors
    Noise(x::Vector{T}, y::Vector{T})  where {T<:AbstractFloat}
    Noise(x::Vector{T}, y::Vector{T}, n::Integer)  where {T<:AbstractFloat}
    Noise(s::Curve{T}) where {T<:AbstractFloat}
    Noise(s::Curve{T}, n::Integer) where {T<:AbstractFloat}

"""
struct Noise{T} <: AbstractCurve{T}
    x::Vector{T}
    y::Vector{T}
    function Noise{T}(x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
        verify_same_length(x, y)
        return new{T}(x, detrend(x, y, 0))
    end
end

Noise(x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat} = Noise{T}(x, y)

function Noise(x::Vector{T}, y::Vector{T}, n::Integer) where {T<:AbstractFloat}
    y = detrend(x, y, n)
    return Noise{T}(x, y)
end

function Noise(s::Curve{T}, n::Integer=0) where {T<:AbstractFloat}
    return Noise(s.x, s.y, n)
end

#######################################################
### Noise models ######################################
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

#######################################################
### Bounds  ###########################################
#######################################################

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
function scale_shift_beta(α::Float64, β::Float64, a::Float64, b::Float64)
    return LocationScale(a, b - a, Beta(α, β))
end

"""
    AbstractUncertainBound

Supertype of uncertain bounds.
"""
abstract type AbstractUncertainBound end

struct LeftRightBound{S1<:ContinuousUnivariateDistribution,
                      S2<:ContinuousUnivariateDistribution} <: AbstractUncertainBound
    left::S1
    right::S2
end

struct WidthBound{T<:ContinuousUnivariateDistribution} <: AbstractUncertainBound
    loc::AbstractFloat
    width::T
end