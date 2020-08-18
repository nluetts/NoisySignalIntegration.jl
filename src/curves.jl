# Defines the Curve and UncertainCurve types for handeling (uncertain) x-y data

"""
    AbstractCurve

Abstract supertype for all curves (xy-data).
Subtypes need to implement fields `x` and `y`.
"""
abstract type AbstractCurve end

Base.eltype(s::AbstractCurve) = eltype(s.x)
Base.length(s::AbstractCurve) = length(s.x)
Base.lastindex(s::AbstractCurve) = length(s)
Base.getindex(s::AbstractCurve, i::Integer) = (s.x[i], s.y[i])
Base.getindex(s::AbstractCurve, r::UnitRange) = typeof(s)(s.x[r], s.y[r])
Base.iterate(s::AbstractCurve, i::Int64) = i > length(s) ? nothing : (s[i], i + 1)
Base.iterate(s::AbstractCurve) = iterate(s, 1)
Base.IteratorSize(itertype::Type{AbstractCurve}) = Base.HasLength()
Base.minimum(s::AbstractCurve) = s[argmin(s.y)]
Base.maximum(s::AbstractCurve) = s[argmax(s.y)]

Base.:(==)(c0::AbstractCurve, c1::AbstractCurve) = c0.x == c1.x && c0.y == c1.y
Base.:(+)(c::AbstractCurve, y) = typeof(c)(c.x, c.y .+ y)
Base.:(-)(c::AbstractCurve, y) = typeof(c)(c.x, c.y .- y)
Base.:(*)(c::AbstractCurve, y) = typeof(c)(c.x, c.y .* y)
Base.:(/)(c::AbstractCurve, y) = typeof(c)(c.x, c.y ./ y)
Base.:(+)(y, c::AbstractCurve) = typeof(c)(c.x, y .+ c.y)
Base.:(-)(y, c::AbstractCurve) = typeof(c)(c.x, y .- c.y)
Base.:(*)(y, c::AbstractCurve) = typeof(c)(c.x, y .* c.y)
Base.:(/)(y, c::AbstractCurve) = typeof(c)(c.x, y ./ c.y)
Base.vcat(c0::AbstractCurve, c1::AbstractCurve) = typeof(c0)(vcat(c0.x, c1.x), vcat(c0.y, c1.y))

function Base.show(io::IO, c::AbstractCurve)
    println(io, "$(typeof(c)), $(length(c)) datapoints")
    if length(c) > 10
        ps = [c[1:5]; c[end-4:end]]
        indicator = "    ⋮\n"
    else
        ps = c
        indicator = ""
    end
    for (i, p) in enumerate(ps)
        println(io, p)
        i == 5 ? print(io, indicator) : nothing
    end
end


"""
    Curve{T} <: AbstractCurve

Datatype holding x-y data.
x and y vectors have to have the same length.

# Fields
- `x::Vector{T}` : x-grid
- `y::Vector{T}` : y-values

# Constructors

    Curve(x::Vector{T}, y::Vector{T})
    Curve(y::Vector{T})

# Examples

```jldoctest
julia> curve = Curve([2.0, 3.0], [4.0, 5.0])
Curve{Float64}, 2 datapoints
(2.0, 4.0)
(3.0, 5.0)


julia> curve = Curve([2.0, 3.0], [5.0])
ERROR: ArgumentError: x and y need to have the same length.
[...]


julia> curve = Curve([4, 5, 6])
Curve{Int64}, 3 datapoints
(1, 4)
(2, 5)
(3, 6)
```
"""
struct Curve{T} <: AbstractCurve
    x::Vector{T}
    y::Vector{T}
    function Curve{T}(x::Vector{T}, y::Vector{T}) where {T}
        verify_same_length(x, y)
        return new{T}(x, y)
    end
end
Curve(x::Vector{T}, y::Vector{T}) where T = Curve{T}(x, y)
Curve(y::Vector{T}) where T = Curve(collect(T, 1:length(y)), y)

"""
    UncertainCurve{T} <: AbstractCurve

Datatype holding x-y data where y-data holds a vector of particles to model uncertainty.
x and y vectors have to have the same length.

# Fields
- `x::Vector{T}` : x-grid
- `y::Vector{Particles{T, N}}` : uncertain y-values

# Constructors

    UncertainCurve(x::Vector{T}, y::Vector{Particles{N, T}})
    UnertainCurve(y::Vector{T})

# Examples

```jldoctest
julia> using MonteCarloMeasurements


julia> uc = UncertainCurve([2.0, 3.0], [4.0 ± 1.0 , 5.0 ± 1.0])
UncertainCurve{Float64,2000}, 2 datapoints
(2.0, 4.0 ± 1.0)
(3.0, 5.0 ± 1.0)


julia> uc = UncertainCurve([2.0], [4.0 ± 1.0 , 5.0 ± 1.0])
ERROR: ArgumentError: x and y need to have the same length.
[...]


julia> uc = UncertainCurve([3.0 ± 1.0, 4.0 ± 1.0 , 5.0 ± 1.0])
UncertainCurve{Float64,2000}, 3 datapoints
(1.0, 3.0 ± 1.0)
(2.0, 4.0 ± 1.0)
(3.0, 5.0 ± 1.0)
```
"""
struct UncertainCurve{T, N} <: AbstractCurve
    x::Vector{T}
    y::Vector{Particles{T, N}}
    function UncertainCurve{T, N}(x::Vector{T}, y::Vector{Particles{T, N}}) where {T, N}
        verify_same_length(x, y)
        return new{T, N}(x, y)
    end
end
UncertainCurve(x::Vector{T}, y::Vector{Particles{T, N}}) where {T, N} = UncertainCurve{T, N}(x, y)
UncertainCurve(y::Vector{Particles{T, N}}) where {T, N} = UncertainCurve(collect(T, 1:length(y)), y)


get_draw(n, p::Particles) = p.particles[n]

"""
    get_draw(n, uc::UncertainCurve)

Retrieve the `n`th sample of the samples stored in `UncertainCurve` `uc`.
"""
get_draw(n, uc::UncertainCurve) = Curve(uc.x, get_draw.(n, uc.y))


Statistics.mean(uc::UncertainCurve) = Curve(uc.x, Statistics.mean(uc.y))


# utility functions

"""
    crop(s::AbstractCurve, left, right)

Crop a slice from a curve and return a new curve.

# Examples


"""
function crop(s::AbstractCurve, left, right)
    T = eltype(s.x)
    S = typeof(s)
    i = searchsortedlast(s.x, T(left))
    j = searchsortedfirst(s.x, T(right))
    return S(s.x[i:j], s.y[i:j])
end