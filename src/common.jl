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
    left = p - w/2
    right = p + w/2
    l = searchsortedfirst(x, left)   # left index:  x[l] <= p - w
    r = searchsortedlast(x, right)  # right index: x[r] >= p + w
    
    m = l + argmax(view(y, l:r)) - 1 # index of local maximum at x[m]
    return [x[m] - w/2, x[m] + w/2]
end


"""Calculate area of single trapezoid."""
@inline singletrapz(x₀, x₁, y₀, y₁) = 0.5 * abs(x₁ - x₀) * (y₁ + y₀)


"""Linearly interpolate y-value at position x between two points (x₀, y₀) and (x₁, y₁)."""
@inline function lininterp(x, x₀, x₁, y₀, y₁)
    δx = x₁ - x₀
    y = (y₁*(x - x₀) + y₀*(x₁ - x)) / δx
    return y
end

"""Linearly interpolate y-value at position `x` for x,y data; `xs` has to be sorted."""
@inline function lininterp(x::T, xs::Vector{T}, ys::Vector{S}) where {T <: Number, S <: Number}
    i = searchsortedfirst(xs, x)
    (i < 2 || i > length(xs)) && throw(error("x is outside the support: $x ∉ ($(minimum(xs)), $(maximum(xs))]."))
    return lininterp(x, xs[i-1], xs[i], ys[i-1], ys[i])
end

"""
    trapz(x::AbstractArray{T}, y::AbstractArray{T}, left, right; subtract_baseline=true) where {T<:AbstractFloat}

Integrate vector `y` in interval [`left`, `right`] using trapezoidal integration.

# Keyword Arguments

`subtract_baseline`: If `true`, subtract a linear baseline defined by data points at `x = left` and `x = right`.

# Notes

`left` and `right` must support conversion to the datatype `T`.

If `left` and `right` do not fall on the `x`-grid, additional data points will be interpolated linearly.
(i.e. the width of the first and last trapezoid will be somewhat smaller).

If `left` and/or `right` falls outside the `x`-range, the integration window will be cropped
to the available range.

# Examples

```jldoctest
julia> x = collect(Float64, 0:10);

julia> y = ones(Float64, 11);


julia> trapz(x, y, 1, 3, subtract_baseline=false)
2.0


julia> trapz(x, y, 1, 3) # subtract_baseline is true by default
0.0


julia> trapz(x, y, -1, 11, subtract_baseline=false) # at most, we can integrate the available x-range, 0 to 10
10.0


julia> trapz(x, y, -10, 20, subtract_baseline=false)
10.0


julia> trapz(x, y, 1.1, 1.3, subtract_baseline=false) ≈ 0.2 # if we integrate "between" the grid, data points are interpolated
true
```
"""
function trapz(x::AbstractArray{T}, y::AbstractArray{T}, left::T, right::T; subtract_baseline=true) where {T<:AbstractFloat}

    left, right = left < right ? (left, right) : (right, left)

    A = zero(eltype(y)) # Area to be returned
    
    N = length(x)
    N != length(y) && throw(ArgumentError("Data arrays `x` and `y` must have the same length."))
    N < 2 && return A
    
    yl = yr = nothing # the y-values of the left and right integration bound, to be interpolated
    lastiter = false
    j = 2

    while j <= N

        x₀ = x[j-1]
        x₁ = x[j]
        y₀ = y[j-1]
        y₁ = y[j]
        
        if x₁ <= left
            j += 1
            continue
        elseif yl === nothing
            # this will only run once, when we enter the integration window
            # test whether x₀ should be replaced by `left`
            if x₀ < left
                y₀ = lininterp(left, x₀, x₁, y₀, y₁)
                x₀ = left
            else
                # this case means that `left` <= x[1]
                left = x₀
            end
            yl = y₀
        end
        
        # test whether x₁ should be replaced by `right`
        if x₁ >= right
            # we move out of the integration window
            yr = x₁ == right ? y₁ : lininterp(right, x₀, x₁, y₀, y₁)
            x₁ = right
            y₁ = yr
            lastiter = true # we shall break the loop after this iteration
        end
        
        A += singletrapz(x₀, x₁, y₀, y₁)
        
        lastiter && break
        
        j += 1
    end

    if subtract_baseline
        if yr === nothing
            # this case means that right > x[end]
            right = x[N]
            yr = y[N]
        end
        return A - singletrapz(left, right, yl, yr)
    else
        return A
    end
end

function trapz(x::AbstractArray{T}, y::AbstractArray{T}, left, right; subtract_baseline=true) where {T<:AbstractFloat}
    left = T(left)
    right = T(right)
    return trapz(x, y, left, right, subtract_baseline=subtract_baseline)
end


"""
    @samples n::Integer e::Expression

Create a `Particles` object with `n` samples using a `±` or `..` expression.

# Examples

```jldoctest
julia> using NoisySignalIntegration

julia> @samples 9999 1 ± 1
Particles{Float64,9999}
 1.0 ± 1.0


julia> @samples 9999 1 .. 2
Particles{Float64,9999}
 1.5 ± 0.289
```
"""
macro samples(n::Integer, e::Expr)
    err = ArgumentError("Expression not understood.")
    if length(e.args) != 3
        return :( throw($err) )
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
        return :( throw($err) )
    end
end