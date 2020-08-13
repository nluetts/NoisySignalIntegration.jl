function get_left_right_points(x::AbstractArray{T}, y::AbstractArray{T}, left::T, right::T; subtract_baseline=true) where {T<:AbstractFloat}
    
    left, right = left < right ? (left, right) : (right, left)

    N = length(x)
    N != length(y) && throw(ArgumentError("Data arrays `x` and `y` must have the same length."))
    N < 2 && throw(ArgumentError("`x`and `y` must each hold at least two elements."))
    
    l =  r = xl = xr = yl = yr = nothing
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
            if x₀ < left
                xl = left
                yl = lininterp(left, x₀, x₁, y₀, y₁)
                l = j - 1
            elseif x₀ == left
                xl = left
                yl = y₀
                l = j
            else
                # this case means that `left` <= x[1]
                xl = x₀
                yl = y₀
                l = 0
            end
        end
        
        if x₁ >= right
            # we move out of the integration window
            yr = x₁ == right ? y₁ : lininterp(right, x₀, x₁, y₀, y₁)
            xr = right
            r = j + 1
            break
        end
        j += 1
    end

    if yr === nothing
        # this case means that right > x[end]
        xr = x[N]
        yr = y[N]
        r = N + 1
    end
    return l, r, xl, xr, yl, yr
end

@recipe function plot_recipe(crv::AbstractCurve)
    xguide --> "x"
    yguide --> "y"
    return crv.x, crv.y
end


# --------------------------------------
# enable plotting of noise sample draws
# --------------------------------------

@recipe function plot_recipe(ns::Curve, uc::UncertainCurve, noise_samples=4)
    noise_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))

    legend := :none
    layout := @layout grid(noise_samples + 1, 1)
    link := :both
    
    delete!(plotattributes, :noise_samples)
    
    for i ∈ 0:noise_samples
        @series begin
            yguide := i == 0 ? "input" : "sample $(i)"
            subplot := i + 1
            if i == 0
                ns.x, ns.y
            else
                uc.x, get_draw(i, uc)
            end
        end
    end
end


@recipe function plot_recipe(crv::Curve{T}, left, right; subtract_baseline=true) where T
    fillrange := 0
    fillalpha --> 0.5
    fillcolor --> :orange
    linewidth --> 0
    label     --> nothing

    left = T(left)
    right = T(right)

    l, r, xl, xr, yl, yr = get_left_right_points(crv.x, crv.y, left, right)
    if subtract_baseline
        x = [xl; crv.x[l+1:r-1]; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; yl]
    else
        x = [xl; crv.x[l+1:r-1]; xr; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; 0;  0 ]
    end
    return x, y
end


# -----------------------------------
# enable plotting of auto covariance
# -----------------------------------
struct AutoCovPlot end # just for dispatch

@recipe function plot_repice(ns::NoiseSample, nm::MvGaussianNoiseModel, ::Type{AutoCovPlot})
    xguide := "lag"
    yguide := "auto-covariance"

    lags, acov = estimate_autocov(ns)

    @series begin
        label := "estimate"
        seriescolor --> :blue
        lags, acov
    end

    @series begin
        label := @sprintf "fit (α = %.3e, λ = %.3e)" nm.α nm.λ
        seriescolor --> :red
        lags, gauss_kernel(lags, [nm.α, nm.λ])
    end
end

plot_autocov(ns::NoiseSample, nm::MvGaussianNoiseModel; kw...) = plot(ns, nm, AutoCovPlot; kw...)