using NoisySignalIntegration
using NoisySignalIntegration: get_draw, _local_baseline, lininterp
using MonteCarloMeasurements
using Plots
using Random
using Debugger

module FWHMmod

struct FWHM{T}
    full_width::T
    _left::T
    _right::T
    _xmax::T
    _half_maximum::T
    _half_maximum_offset::T
end

end

function mc_fwhm(uc::UncertainCurve{T,N}, bnds::Vector{UncertainBound{T,M}}; local_baseline=false) where {T,M,N}

    M != N && error("Samples sizes incompatible")

    widths = Array{FWHMmod.FWHM{T}}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Processing draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            xₗ, xᵣ = get_draw(i, b)
            if local_baseline
                baseline = _local_baseline(cᵢ.x, cᵢ.y, xₗ, xᵣ, b)
                widths[i, j] = fwhm(cᵢ, xₗ, xᵣ, baseline)
            else
                widths[i, j] = fwhm(cᵢ, xₗ, xᵣ)

            end
        end
    end
    return widths
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
    xₗ, xᵣ, yₗ, yᵣ = isnothing(baseline) ? (left, right, zero(T), zero(T)) : baseline

    # index of peak maximum
    imx = argmax(ys)
    xmx = xs[imx]
    # y-value of baseline at maximum
    ybx = lininterp(xmx, xₗ, xᵣ, yₗ, yᵣ)
    # half-maximum
    hm = (maximum(ys) - ybx)/2

    left = NaN
    right = NaN
    # find left position
    for i in 1:(imx - 1)
        yi = ys[i] - yₗ
        yj = ys[i + 1] - yₗ
        @bp
        if yj >= hm && yi < hm 
            # find x-value that corresponds best to max/2
            left = lininterp(hm, yi, yj, xs[i], xs[i + 1])
            break
        end
    end
    # find right position
    for i in imx:(length(ys) - 1)
        yi = ys[i] - yᵣ
        yj = ys[i + 1] - yᵣ
        if yi >= hm && yj < hm
            # find x-value that corresponds best to max/2
            right = lininterp(hm, yi, yj, xs[i], xs[i + 1])
            break
        end
    end

    full_width = right - left
    # for plotting purposes, we return bounds and half-maximum as well as fwhm
    FWHMmod.FWHM(full_width, left, right, xmx, hm, ybx)
end

function plot_fwhm!(p, f::FWHMmod.FWHM{T}) where {T<:Number}
    plot!(p, [f._xmax, f._xmax], [f._half_maximum_offset, f._half_maximum_offset + f._half_maximum], color=:gray, label=nothing, linestyle=:dot)
    plot!(p, [f._xmax, f._xmax], [f._half_maximum_offset + f._half_maximum, f._half_maximum_offset + f._half_maximum * 2], color=:gray, label=nothing, linestyle=:dot)
    plot!(p, [f._left, f._right], [f._half_maximum + f._half_maximum_offset, f._half_maximum + f._half_maximum_offset], color=:gray, label=nothing)
end


Random.seed!(42)
spectrum = NoisySignalIntegration.testdata_1()
slice_bands = crop(spectrum,  5.0,  40.0)
slice_noise = crop(spectrum, 40.0, 100.0)

noise = NoiseSample(slice_noise, 3)
nm = fit_noise(noise)
uncertain_spectrum = add_noise(slice_bands, nm)
position = [15.0, 30.0]
# widths will fall in the range 2 to 3, with a maximum at 2.5
width_distribution = scale_shift_beta(2, 2, 3, 4)
# define a "width bound"
bds = UncertainBound(position, width_distribution, uncertain_spectrum)
#@run mc_fwhm(uncertain_spectrum, bds, local_baseline=true)
p = plot(spectrum)
xlims!(0, 50)

xl, xr = get_draw(1, bds[2])
baseline = _local_baseline(spectrum.x, spectrum.y, xl, xr, bds[2])
width = fwhm(spectrum, xl, xr, baseline)
plot_fwhm!(p, width)
xl, xr = get_draw(1, bds[1])
baseline = _local_baseline(spectrum.x, spectrum.y, xl, xr, bds[2])
width = fwhm(spectrum, xl, xr, baseline)
plot_fwhm!(p, width)

