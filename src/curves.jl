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

stitch(curves...) = reduce(vcat, curves)

function Base.sort(c::T) where {T <: AbstractCurve}
    x, y = c.x, c.y
    ind = sortperm(x)
    return T(x[ind], y[ind])
end

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
- `x :: Vector{T}` : x-grid
- `y :: Vector{T}` : y-values


# Constructors

    Curve(x, y)
    Curve(y)

# Note

Both `x` and `y` must be convertabile to vectors with the same element type.

# Examples

```jldoctest
julia> Curve([2.0, 3.0], [4.0, 5.0])
Curve{Float64}, 2 datapoints
(2.0, 4.0)
(3.0, 5.0)


julia> Curve([2.0, 3.0], [5.0])
ERROR: ArgumentError: x and y need to have the same length.
[...]


julia> Curve([4, 5, 6])
Curve{Int64}, 3 datapoints
(1, 4)
(2, 5)
(3, 6)


julia> Curve(1:3, [1., 2., 3.])
Curve{Float64}, 3 datapoints
(1.0, 1.0)
(2.0, 2.0)
(3.0, 3.0)
```
"""
struct Curve{T} <: AbstractCurve
    x::Vector{T}
    y::Vector{T}
    function Curve{T}(x, y) where {T <: Number}
        verify_same_length(x, y)
        return new{T}(x, y)
    end
end

function Curve(y)
    T = eltype(y)
    return Curve{T}(collect(T, 1:length(y)), y)
end

Curve(x::Vector{T}, y::Vector{T}) where T = Curve{T}(x, y)

function Curve(x, y)
    verify_same_length(x, y)
    xypairs = promote.(x, y)
    x = map(first, xypairs)
    y = map(last, xypairs)
    T = eltype(y)
    Curve{T}(x, y)
end

"""
    UncertainCurve{T} <: AbstractCurve

Datatype holding x-y data where y-data holds a vector of particles to model uncertainty.
`x` and `y` vectors have to have the same length.

# Fields
- `x::Vector{T}` : x-grid
- `y::Vector{Particles{T, N}}` : uncertain y-values

# Constructors

    UncertainCurve(x::Vector{T}, y::Vector{Particles{N, T}})
    UnertainCurve(y::Vector{T})

# Examples

```jldoctest
julia> using MonteCarloMeasurements


julia> UncertainCurve([2.0, 3.0], [4.0 ± 1.0 , 5.0 ± 1.0])
UncertainCurve{Float64, 2000}, 2 datapoints
(2.0, 4.0 ± 1.0)
(3.0, 5.0 ± 1.0)


julia> UncertainCurve([2.0], [4.0 ± 1.0 , 5.0 ± 1.0])
ERROR: ArgumentError: x and y need to have the same length.
[...]


julia> UncertainCurve([3.0 ± 1.0, 4.0 ± 1.0 , 5.0 ± 1.0])
UncertainCurve{Float64, 2000}, 3 datapoints
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


"""
    get_draw(n, p::Particles)

Retrieve the `n`th particle of the particles stored in `p`.
"""
get_draw(n, p::Particles) = p.particles[n]

"""
    get_draw(n, uc::UncertainCurve)

Retrieve the `n`th sample of the samples stored in `UncertainCurve` `uc`.
"""
get_draw(n, uc::UncertainCurve) = Curve(uc.x, get_draw.(n, uc.y))

"""
    mean(uc::UncertainCurve)

Retrieve the mean of the `UncertainCurve` `uc`.
"""
Statistics.mean(uc::UncertainCurve) = Curve(uc.x, Statistics.mean(uc.y))


MonteCarloMeasurements.:(±)(c::Curve, v) = UncertainCurve(c.x, c.y ± v)

# utility functions

"""
    crop(s::AbstractCurve, left, right)

Crop a slice from a curve between `left` and `right` and return a new curve.

# Examples

```jldoctest
julia> c = Curve(1:0.5:10, 11:0.5:20)
Curve{Float64}, 19 datapoints
(1.0, 11.0)
(1.5, 11.5)
(2.0, 12.0)
(2.5, 12.5)
(3.0, 13.0)
    ⋮
(8.0, 18.0)
(8.5, 18.5)
(9.0, 19.0)
(9.5, 19.5)
(10.0, 20.0)


julia> crop(c, 2, 3)
Curve{Float64}, 3 datapoints
(2.0, 12.0)
(2.5, 12.5)
(3.0, 13.0)


julia> crop(c, 2.1, 3.5)
Curve{Float64}, 3 datapoints
(2.5, 12.5)
(3.0, 13.0)
(3.5, 13.5)


julia> crop(c, 2.1, 3.9)
Curve{Float64}, 3 datapoints
(2.5, 12.5)
(3.0, 13.0)
(3.5, 13.5)


julia> crop(c, 2.1, 4.0)
Curve{Float64}, 4 datapoints
(2.5, 12.5)
(3.0, 13.0)
(3.5, 13.5)
(4.0, 14.0)
```
"""
function crop(s::AbstractCurve, left, right)
    T = eltype(s.x)
    S = typeof(s)
    i = searchsorted(s.x, T(left))  |> first
    j = searchsorted(s.x, T(right)) |> last
    return S(s.x[i:j], s.y[i:j])
end