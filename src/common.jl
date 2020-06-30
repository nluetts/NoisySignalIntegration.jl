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
    r = searchsortedfirst(x, right)  # right index: x[r] >= p + w
    
    m = l + argmax(view(y, l:r)) - 1 # index of local maximum at x[m]
    return [x[m] - w/2, x[m] + w/2]
end

function get_linear_param(x₀, Δx, y₀, y₁)
    m = (y₁-y₀)/Δx
    a = y₀-x₀*m
    a, m # offset and slope
end

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

    return δx, l, r, x_ll, x_lr, x_rl, x_rr, y_ll, y_lr, y_rl, y_rr, x_l, x_r, y_l, y_r
end

"""
    function trapz(x, y, left, right))

Integrate vector `x` in interval [`left`, `right`] using trapezoidal integration
after subtracting a baseline defined by data points at `x = left, right`.
"""
function trapz(x, y, left, right)

    δx, l, r, x_ll, x_lr, x_rl, x_rr, y_ll, y_lr, y_rl, y_rr, x_l, x_r, y_l, y_r = get_integration_bounds(x, y, left, right)

    # integral correction for interpolated bounds
    A = 0.5 * ((x_lr - x_l)*(y_lr + y_l)
             + (x_r - x_rl)*(y_rl + y_r));

    # trapezoidal integration
    y₀ = y_lr
    for j = (l+2):(r-1)
        y₁ = y[j]
        A += 0.5*δx*(y₀ + y₁)
        y₀ = y₁
    end
    A - 0.5*(x_r-x_l)*(y_l+y_r)
end;