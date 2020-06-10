"""
    function get_linear_param(x₀, Δx, y₀, y₁)
"""
function get_linear_param(x₀, Δx, y₀, y₁)
    m = (y₁-y₀)/Δx
    a = y₀-x₀*m
    a, m # offset and slope
end;


"""
    function trapz(x, y, left, right))

Integrate vector `x` in interval [`left`, `right`] using trapezoidal integration
after subtracting a baseline defined by data points at `x = left, right`.
"""
function trapz(x, y, left, right)

    δx = x[2]-x[1] # x increment (x must be evenly spaced!)

    # find indices of interval [left, right]
    l = searchsortedlast(x, left)   # left index:  x[l] <= left
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