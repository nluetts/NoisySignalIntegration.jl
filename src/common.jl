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
    detrend(x, y, poly_order)

Subtract polynomial from y data.
"""
detrend(x, y, poly_order) = y - fit(x, y, poly_order).(x)

"""
    left_right_from_peak(x, y, p, w)

Find peak in interval `p ± w/2` and return symmetric bounds
left and right from peak with width `w`.
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

"""
    function trapz(x, y, left, right))

Integrate vector `x` in interval [`left`, `right`] using trapezoidal integration
after subtracting a baseline defined by data points at `x = left, right`.
"""
@inline singletrapz(x₀, x₁, y₀, y₁) = 0.5 * abs(x₁ - x₀) * (y₁ + y₀)
@inline function lininterp(x, x₀, x₁, y₀, y₁)
    δx = x₁ - x₀
    y = (y₁*(x - x₀) + y₀*(x₁ - x)) / δx
    return y
end


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