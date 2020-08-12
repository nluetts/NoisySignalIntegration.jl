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

struct UncertainBound{T, N}
    left::Particles{T, N}
    right::Particles{T, N}
end

get_draw(n, p::Particles) = p.particles[n]
get_draw(n, uc::UncertainCurve) = [get_draw(n, yᵢ) for yᵢ in uc.y]
get_draw(n, ub::UncertainBound) = [get_draw(n, ub.left), get_draw(n, ub.right)]

# create a left/right bound
function UncertainBound(left::S, right::T, N::Int=10_000) where {S <: ContinuousUnivariateDistribution, T <: ContinuousUnivariateDistribution}
    left  = Particles(N, left)
    right = Particles(N, right)
    return UncertainBound(left, right)
end

# create a width bounds with correlated widths
function UncertainBound(
    pos::Vector{T},
    width::ContinuousUnivariateDistribution,
    uc::UncertainCurve{T, N}
) where {T, N}

    M = length(pos)
    left = Array{T}(undef, N, M)
    right = Array{T}(undef, N, M)
    width = Particles(N, width)

    for i ∈ 1:M
        pᵢ = pos[i]
        for j ∈ 1:N
            cⱼ = get_draw(j, uc) # sample j from curve
            wⱼ = get_draw(j, width) # sample j from width
            left[i, j], right[i, j] = left_right_from_peak(uc.x, cⱼ, pᵢ, wⱼ)
        end
    end
    
    return [UncertainBound(Particles(left[i, :]), Particles(right[i, :])) for i in 1:M]
end

# create single width bound
function UncertainBound(
    pos::T,
    width::ContinuousUnivariateDistribution,
    uc::UncertainCurve{T, N}
) where {T, N}
    bnd = UncertainBound([pos], width, uc)
    return bnd[1]
end

ScaledShiftedBeta = LocationScale{Float64, Beta{Float64}}

"""
    scale_shift_beta(α, β, a, b)

Create a scaled and shifted `Beta(α, β)` distribution.
Samples fall in the interval [`a`, `b`].
"""
function scale_shift_beta(α, β, a, b)
    return LocationScale(a, b - a, Beta(α, β))
end