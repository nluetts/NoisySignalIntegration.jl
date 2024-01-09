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
Do nothing if x and y have same length, otherwise throw ArgumentError.
"""
function verify_same_length(x::AbstractArray, y::AbstractArray)
    length(x) == length(y) || throw(ArgumentError("x and y need to have the same length."))
    return nothing
end


"""
Find peak in interval `p ± w/2` and return (speak position - `w`/2, peak position + `w`/2).
"""
function left_right_from_peak(x, y, p, w)
    # find indices of interval [left, right]
    left = p - w / 2
    right = p + w / 2
    l = searchsortedfirst(x, left)   # left index:  x[l] <= p - w
    r = searchsortedlast(x, right)  # right index: x[r] >= p + w

    m = l + argmax(view(y, l:r)) - 1 # index of local maximum at x[m]
    return [x[m] - w / 2, x[m] + w / 2]
end


"""Calculate area of single trapezoid."""
@inline singletrapz(x₀, x₁, y₀, y₁) = 0.5 * abs(x₁ - x₀) * (y₁ + y₀)


"""Linearly interpolate y-value at position x between two points (x₀, y₀) and (x₁, y₁)."""
@inline function lininterp(x, x₀, x₁, y₀, y₁)
    δx = x₁ - x₀
    y = (y₁ * (x - x₀) + y₀ * (x₁ - x)) / δx
    return y
end

"""Linearly interpolate y-value at position `x` for x,y data; `xs` has to be sorted."""
@inline function lininterp(x::T, xs::Vector{T}, ys::Vector{S}) where {T<:Number,S<:Number}
    i = searchsortedfirst(xs, x)
    (i < 2 || i > length(xs)) && throw(error("x is outside the support: $x ∉ ($(minimum(xs)), $(maximum(xs))]."))
    return lininterp(x, xs[i-1], xs[i], ys[i-1], ys[i])
end


"""
    @samples n::Integer e::Expression

Create a `Particles` object with `n` samples using a `±` or `..` expression.

# Examples

```jldoctest
julia> using NoisySignalIntegration

julia> @samples 9999 1 ± 1
MonteCarloMeasurements.Particles{Float64, 9999}
 1.0 ± 1.0


julia> @samples 9999 1 .. 2
MonteCarloMeasurements.Particles{Float64, 9999}
 1.5 ± 0.289
```
"""
macro samples(n::Integer, e::Expr)
    err = ArgumentError("Expression not understood.")
    if length(e.args) != 3
        return :(throw($err))
    end

    op, x, y = e.args

    if op == :±
        return quote
            a = $(esc(x))
            b = $(esc(y))
            Particles($n, Normal(a, b))
        end
    elseif op == :..
        return quote
            a = $(esc(x))
            b = $(esc(y))
            Particles($n, Uniform(a, b))
        end
    else
        return :(throw($err))
    end
end
