function baseline_from_points(x1, x2, y1, y2)
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    @assert b ≈ y2 - m * x2 rtol=1e-10
    m, b
end


function band_center(
    curve::Curve{T},
    left::T,
    right::T,
    baseline::Union{Nothing, Tuple{T, T}}=nothing
) where {T<:AbstractFloat}
    sum_intensity = zero(T)
    sum_weighted_position = zero(T)
    m, b = isnothing(baseline) ? (zero(T), zero(T)) : baseline
    for (xi, yi) in (curve, left, right)
        sum_intensity += yi - m * xi - b
        sum_weighted_position += xi * (yi - m * xi - b)
    end
    sum_weighted_position / sum_intensity
end

# this will behave like `subtract_baseline`, i.e. spanning a linear
# baseline from the start to the endpoint
function band_center(
    curve::Curve{T},
    left::T,
    right::T,
    baseline::Bool
) where {T<:AbstractFloat}
    yleft = lininterp(left, curve)
    yright = lininterp(right, curve)
    if baseline
        m, b = baseline_from_points(left, right, yleft, yright)
        band_center(curve, left, right, (m, b))
    else
        band_center(curve, left, right, nothing)
    end
end


"""
    mc_bandcenter(uc::UncertainCurve{T, N}, bnds::Vector{UncertainBound{T, M}})
    mc_bandcenter(uc::UncertainCurve{T, N}, bnds::UncertainBound{T, M})

Determines band centers from uncertain curve using uncertain bound(s).
Band centers are calculated as the position weighted by intensity (Σᵢxᵢyᵢ)/(Σᵢyᵢ).
Much like with the function `trapz` when integrating, data points are interpolated linearly if
the integration bounds happen to fall between the grid spanned by the x-values.

Returns one or an array of band centers of type `Particles{T, N}`.

# Keyword arguments

`subtract_baseline` (deprecated in favor of `local_baseline`):
If true, for each draw a local linear baseline defined by the integration window start and end point will be subtracted.

`local_baseline`:
If true, for each draw a local linear baseline defined by the integration window start and end point will be subtracted.
The y-values of the start and end point are derived from a weighted average over the start and end point distributions, see 
[the documentation](https://nluetts.github.io/NoisySignalIntegration.jl/dev/baseline/#Build-in) for further information.
"""
function mc_bandcenter(uc::UncertainCurve{T,N}, bnds::Vector{UncertainBound{T,M}}; subtract_baseline=false, local_baseline=false) where {T,M,N}

    M != N && error("Samples sizes incompatible")
    subtract_baseline && @warn("subtract_baseline keyword argument is deprecated, use local_baseline instead.")
    (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

    centers = Array{T}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Processing draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            xₗ, xᵣ = get_draw(i, b)
            if local_baseline
                baseline = baseline_from_points(_local_baseline(cᵢ.x, cᵢ.y, xₗ, xᵣ, b)...)
                centers[i, j] = band_center(cᵢ, xₗ, xᵣ, baseline)
            elseif subtract_baseline
                centers[i, j] = band_center(cᵢ, xₗ, xᵣ, true)
            else
                centers[i, j] = band_center(cᵢ, xₗ, xᵣ)

            end
        end
    end
    return [Particles(centers[:, i]) for i in 1:size(centers)[2]]
end

function mc_bandcenter(uc::S, bnd::T; subtract_baseline=false, local_baseline=false) where {S<:UncertainCurve,T<:UncertainBound}
    return mc_bandcenter(uc, [bnd]; subtract_baseline=subtract_baseline, local_baseline=local_baseline)[1]
end
