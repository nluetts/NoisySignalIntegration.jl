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

function get_linear_param(x₀, Δx, y₀, y₁)
    m = (y₁-y₀)/Δx
    a = y₀-x₀*m
    return a, m # offset and slope
end


# this horrible function only exists because the bounds are needed
# for plotting and wrapping them in their own object caused a big
# decrease in performance 
@inline function get_integration_bounds(x, y, left, right)
    δx = x[2]-x[1] # x increment (x must be evenly spaced!)

    # find indices of interval [left, right]
    l = searchsortedlast(x, left)  # left index:  x[l] <= left
    r = searchsortedfirst(x, right) # right index: x[r] >= right

    # boundary data points
    x_ll = x[l]
    x_lr = x[l+1]
    x_rl = x[r-1]
    x_rr = x[r]
    y_ll = y[l]
    y_lr = y[l+1]
    y_rl = y[r-1]
    y_rr = y[r]

    # linearly interpolated boundary data points
    x_l = left
    x_r = right
    y_l = begin
        a_l, m_l = get_linear_param(x_ll, δx, y_ll, y_lr)
        a_l + x_l*m_l
    end
    y_r = begin
        a_r, m_r = get_linear_param(x_rl, δx, y_rl, y_rr)
        a_r + x_r*m_r
    end

    # yes, one should not return a 15-tuple
    return δx, l, r, x_ll, x_lr, x_rl, x_rr, y_ll, y_lr, y_rl, y_rr, x_l, x_r, y_l, y_r
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

function trapz(x::AbstractArray{T}, y::AbstractArray{T}, left::T, right::T) where {T<:AbstractFloat}

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

    # subtract baseline
    if yr === nothing
        # this case means that right > x[end]
        right = x[N]
        yr = y[N]
    end

    A_baseline = singletrapz(left, right, yl, yr)

    return A - A_baseline
end

function trapz(x::AbstractArray{T}, y::AbstractArray{T}, left, right) where {T<:AbstractFloat}
    left = T(left)
    right = T(right)
    return trapz(x, y, left, right)
end