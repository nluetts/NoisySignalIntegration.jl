# Defines the Curve and UncertainCurve types for handeling (uncertain) x-y data

"""
    AbstractCurve

Abstract supertype for all curves (xy-data).
Subtypes need to implement fields `x` and `y`.
"""
abstract type AbstractCurve end

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
    println("$(typeof(c)), $(length(c)) datapoints")
    if length(c) > 10
        ps = [c[1:5]; c[end-4:end]]
        indicator = "    ⋮\n"
    else
        ps = c
        indicator = ""
    end
    for (i, p) in enumerate(ps)
        println(p)#@printf "    %-5e %-5e\n" p...
        i == 5 ? print(indicator) : nothing
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
Curve(y::T) where T = Curve(collect(T, 1:length(y)), y)


struct UncertainCurve{T, N} <: AbstractCurve
    x::Vector{T}
    y::Vector{Particles{T, N}}
end


get_draw(n, p::Particles) = p.particles[n]
get_draw(n, uc::UncertainCurve) = [get_draw(n, yᵢ) for yᵢ in uc.y]


"""
    crop(s::AbstractCurve, left, right)

Crop a slice from a curve and return a new curve.
"""
function crop(s::AbstractCurve, left, right)
    T = eltype(s.x)
    S = typeof(s)
    i = searchsortedlast(s.x, T(left))
    j = searchsortedfirst(s.x, T(right))
    return S(s.x[i:j], s.y[i:j])
end


"""
    detrend(x, y, poly_order)

Subtract polynomial from y data.
"""
detrend(x, y, poly_order) = y - fit(x, y, poly_order).(x)
detrend(c::Curve, poly_order) = Curve(c.x, detrend(c.x, c.y, poly_order))