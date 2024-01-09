function _local_baseline(xs, ys, xₗ, xᵣ, b::UncertainBound)
    ŷₗ = mean(lininterp(q, xs, ys) for q in b._left_quantiles)
    ŷᵣ = mean(lininterp(q, xs, ys) for q in b._right_quantiles)
    x̂ₗ = b._left_quantiles[IDX_BOUND_MEDIAN]
    x̂ᵣ = b._right_quantiles[IDX_BOUND_MEDIAN]
    yₗ = lininterp(xₗ, x̂ₗ, x̂ᵣ, ŷₗ, ŷᵣ)
    yᵣ = lininterp(xᵣ, x̂ₗ, x̂ᵣ, ŷₗ, ŷᵣ)
    return xₗ, xᵣ, yₗ, yᵣ
end


_endpoint_to_endpoint_baseline(xs, ys, xₗ, xᵣ) = (xₗ, xᵣ, lininterp(xₗ, xs, ys), lininterp(xᵣ, xs, ys))

function trapz(crv::Curve{T}, left::T, right::T) where {T<:AbstractFloat}
    left, right = min(left, right), max(left, right)
    area = zero(T)
    for ((xi, yi), (xj, yj)) in zip((crv, left, right), Iterators.drop((crv, left, right), 1))
        area += singletrapz(xi, xj, yi, yj)
    end
    return area
end

"""
    trapz(x::AbstractArray{T}, y::AbstractArray{T}, left, right) where {T<:AbstractFloat}

Integrate vector `y` in interval [`left`, `right`] using trapezoidal integration.

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


julia> trapz(x, y, 1, 3)
2.0


julia> trapz(x, y, -1, 11) # at most, we can integrate the available x-range, 0 to 10
10.0


julia> trapz(x, y, -10, 20)
10.0


julia> trapz(x, y, 1.1, 1.3) ≈ 0.2 # if we integrate "between" the grid, data points are interpolated
true
```
"""
function trapz(xs::AbstractArray{T}, ys::AbstractArray{T}, left::T, right::T) where {T<:AbstractFloat}
    trapz(Curve(xs, ys), left, right)
end

function trapz(x::AbstractArray{T}, y::AbstractArray{T}, left, right) where {T<:AbstractFloat}
    left = T(left)
    right = T(right)
    return trapz(x, y, left, right)
end


"""
    mc_integrate(uc::UncertainCurve{T, N}, bnds::Vector{UncertainBound{T, M}}; intfun=trapz)
    mc_integrate(uc::UncertainCurve{T, N}, bnds::UncertainBound{T, M}; intfun=trapz)

Integrate uncertain curve using uncertain bound(s).

Applies Monte-Carlo integration algorithm to samples stored in `uc` and returns one or an array
of uncertain areas of type `Particles{T, N}`.

# Keyword arguments

`intfun`:
The core integration function that is used to numerically integrate each draw.
Defaults to `NoisySignalIntegration.trapz`.
The function that is used to substitute [`trapz`](@ref) must share its call signature.

`subtract_baseline` (deprecated in favor of `local_baseline`):
If true, for each draw a local linear baseline defined by the integration window start and end point will be subtracted.

`local_baseline`:
If true, for each draw a local linear baseline defined by the integration window start and end point will be subtracted.
The y-values of the start and end point are derived from a weighted average over the start and end point distributions, see 
[the documentation](https://nluetts.github.io/NoisySignalIntegration.jl/dev/baseline/#Build-in) for further information.
"""
function mc_integrate(uc::UncertainCurve{T,N}, bnds::Vector{UncertainBound{T,M}}; intfun=trapz, subtract_baseline=false, local_baseline=false) where {T,M,N}

    M != N && error("Samples sizes incompatible")
    subtract_baseline && @warn("subtract_baseline keyword argument is deprecated, use local_baseline instead.")
    (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

    areas = Array{T}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Integrating draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            xₗ, xᵣ = get_draw(i, b)
            xs, ys = uc.x, cᵢ.y
            areas[i, j] = intfun(xs, ys, xₗ, xᵣ)
            if local_baseline
                areas[i, j] -= singletrapz(_local_baseline(xs, ys, xₗ, xᵣ, b)...)
            end
            if subtract_baseline
                areas[i, j] -= singletrapz(_endpoint_to_endpoint_baseline(xs, ys, xₗ, xᵣ)...)
            end
        end
    end
    return [Particles(areas[:, i]) for i in 1:size(areas)[2]]
end

function mc_integrate(uc::S, bnd::T; intfun=trapz, subtract_baseline=false, local_baseline=false) where {S<:UncertainCurve,T<:UncertainBound}
    return mc_integrate(uc, [bnd]; intfun=intfun, subtract_baseline=subtract_baseline, local_baseline=local_baseline)[1]
end
