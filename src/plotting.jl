# --------------------------------------------
# enable plotting of custom types via recipes
# --------------------------------------------


# --------------------------------------------
# enable plotting of curves and noise samples
# --------------------------------------------
@recipe function plot_recipe(crv::AbstractCurve)
    xguide --> "x"
    yguide --> "y"
    return crv.x, crv.y
end


# --------------------------------------
# enable plotting of noise sample draws
# --------------------------------------

@recipe function plot_recipe(nm::AbstractNoiseModel; grid_points=1000, noise_samples=3)
    legend --> :outertopright
    noise_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    S = sample(nm, grid_points, noise_samples)
    delete!(plotattributes, :grid_points)
    delete!(plotattributes, :noise_samples)
    span = (maximum(S) - minimum(S)) * 1.1
    for (i, s) in enumerate(eachcol(S))
        @series begin
            label := "sample $(i)"
            s .+ (i*span - 1)
        end
    end
end

@recipe function plot_recipe(ns::Noise, nm::AbstractNoiseModel; noise_samples=3)
    noise_samples < 0 && throw(ArgumentError("Number of samples must be > 0."))
    S = sample(nm, length(ns), noise_samples)
    span = (maximum(S) - minimum(S)) * 1.1
    for (i, s) in enumerate(eachcol(S))
        @series begin
            label := "sample $(i)"
            ns.x, s .+ i*span
        end
    end
    seriescolor --> :black
    label := "experimental"
    ns.x, ns.y
end

# -----------------------------------
# enable plotting of auto covariance
# -----------------------------------
struct AutoCovPlot end # just for dispatch

@recipe function plot_repice(ns::Noise, nm::MvGaussianNoiseModel, ::Type{AutoCovPlot})
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

plot_autocov(ns::Noise, nm::MvGaussianNoiseModel; kw...) = plot(ns, nm, AutoCovPlot; kw...)

# --------------------------------------
# enable plotting of integration bounds
# --------------------------------------

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


@recipe function plot_recipe(crv::Curve, left::Float64, right::Float64, blp::BaselinePolicy)
    fillrange := 0
    fillalpha --> 0.5
    fillcolor --> :orange
    linewidth --> 0
    label     --> nothing

    l, r, xl, xr, yl, yr = get_left_right_points(crv.x, crv.y, left, right)
    if blp == SUBTRACT_LOCAL
        x = [xl; crv.x[l+1:r-1]; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; yl]
    else
        x = [xl; crv.x[l+1:r-1]; xr; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; 0;  0 ]
    end
    return x, y
end

@recipe function plot_recipe(crv::Curve, bnds::Vector{T}) where {T <: AbstractUncertainBound}
    for b in bnds
        @series begin
            l, r = typeof(b) <: WidthBoundUnion ? sample(b, crv) : sample(b)
            crv, l, r
        end
    end
    return crv
end

@recipe function plot_recipe(crv::Curve, bnds::Vector{T}, nm::AbstractNoiseModel; samples=3) where {T <: AbstractUncertainBound}
    
    legend --> :best
    background_color_legend --> nothing
    foreground_color_legend --> nothing
    layout := @layout grid(samples+1, 1)
    link := :both

    _, y_min = minimum(crv)
    _, y_max = maximum(crv)

    # plot experimental spectrum with mean bounds
    @series begin
        label := "original data, mean bounds"
        showaxis := :y
        subplot := 1
        crv
    end
    for (i, b) in enumerate(bnds)
        @series begin
            subplot := 1
            if typeof(b) <: WidthBound
                l, r, = left_right_from_peak(crv.x, crv.y, b.loc, mean(b.width))
            elseif typeof(b) == WidthBoundClone
                :fillcolor --> :red
                l, r, = left_right_from_peak(crv.x, crv.y, b.loc, mean(b.reference.width))
            else
                l, r = map(mean, [b.left, b.right])
            end
            crv, l, r, get_baseline_policy(b)
        end
    end

    # plot samples
    for i in 1:samples
        crv_ = crv + (sample(nm, length(crv)))
        @series begin
            label := "sample $i"
            subplot := i + 1
            crv_
        end
    
        for b in bnds
            @series begin
                label := nothing
                subplot := i + 1
                if typeof(b) == WidthBoundClone
                    :fillcolor --> :red
                end
                l, r = typeof(b) <: WidthBoundUnion ? sample(b, crv_) : sample(b)
                crv_, l, r, get_baseline_policy(b)
            end
        end
    end
end