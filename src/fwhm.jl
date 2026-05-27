
"""
This struct is only used internally, for dispatch when plotting.
"""
struct FWHM{T}
    full_width::T
    _local_baseline::Union{Nothing, Tuple{T, T, T, T}}
    _peak_position::T
    _start_position::T
    _half_maximum::T
    _half_maximum_offset::T
end

function value(fwhm::FWHM{T}) :: T where {T}
    fwhm.full_width
end


function mc_fwhm(
    uc::UncertainCurve{T,N},
    bnds::Vector{UncertainBound{T,M}};
    local_baseline=false
) :: Vector{Particles{T}} where {T,M,N}

    M != N && error("Samples sizes of bounds and uncertain curve incompatible ($N != $M)")

    widths = Array{Float64}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Processing draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            xₗ, xᵣ = get_draw(i, b)
            if local_baseline
                baseline = _local_baseline(cᵢ.x, cᵢ.y, xₗ, xᵣ, b)
                widths[i, j] = fwhm(cᵢ, xₗ, xᵣ, baseline) |> value
            else
                widths[i, j] = fwhm(cᵢ, xₗ, xᵣ) |> value

            end
        end
    end

    return [Particles(widths[:, i]) for (i, _ws) in enumerate(eachcol(widths))]
end

function fwhm(
    curve::Curve{T},
    left::T,
    right::T,
    baseline::Union{Nothing, Tuple{T, T, T, T}}=nothing
) where {T<:AbstractFloat}
    mask = curve.x .> left .&& curve.x .<= right
    xs = curve.x[mask]
    ys = curve.y[mask]
    xl, xr, yl, yr = isnothing(baseline) ? (left, right, zero(T), zero(T)) : baseline

    # index and position of peak maximum
    imx = argmax(ys)
    peak_position = xs[imx]
    # y-value of baseline at maximum
    ybx = lininterp(peak_position, xl, xr, yl, yr)
    # half-maximum
    hm = (maximum(ys) - ybx)/2

    x_fwhm_left = NaN
    x_fwhm_right = NaN
    # find left position
    for i in 1:(imx - 1)
        yi = ys[i] - ybx
        yj = ys[i + 1] - ybx
        if yi < hm <= yj 
            # find x-value that corresponds best to max/2
            x_fwhm_left = lininterp(hm, yi, yj, xs[i], xs[i + 1])
            break
        end
    end
    # find right position
    for i in imx:(length(ys) - 1)
        yi = ys[i] - ybx
        yj = ys[i + 1] - ybx
        if yj < hm <= yi 
            # find x-value that corresponds best to max/2
            x_fwhm_right = lininterp(hm, yi, yj, xs[i], xs[i + 1])
            break
        end
    end

    full_width = x_fwhm_right - x_fwhm_left
    FWHM(full_width, baseline, peak_position, x_fwhm_left, hm, ybx)
end
