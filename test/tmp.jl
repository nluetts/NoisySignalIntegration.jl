using NoisySignalIntegration
using NoisySignalIntegration: get_draw, _local_baseline, lininterp, get_left_right_points, _endpoint_to_endpoint_baseline
using MonteCarloMeasurements
using Plots
using RecipesBase
using Random
using Debugger

module FWHMmod

using NoisySignalIntegration: Curve, UncertainCurve, UncertainBound

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

end

using .FWHMmod: FWHM, value

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
    anim = Animation()
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
    FWHMmod.FWHM(full_width, baseline, peak_position, x_fwhm_left, hm, ybx)
end

@recipe function plot_recipe(f::FWHMmod.FWHM{T}) where {T<:Number}
    # vertical line at peak maximum
    @series begin
        color := :gray
        label := nothing
        linestyle := :dot
        [f._peak_position, f._peak_position], [f._half_maximum_offset, f._half_maximum_offset + 2*f._half_maximum]
    end
    # horizontal line marking fwhm
    @series begin
        color := :gray
        label := nothing
        linewidth --> 2.0
        [f._start_position, f._start_position + f.full_width], [f._half_maximum + f._half_maximum_offset, f._half_maximum + f._half_maximum_offset]
    end
    # local baseline, if set
    if !isnothing(f._local_baseline)
        @series begin
            color := :gray
            label := nothing
            linestyle := :dot
            [f._local_baseline[1], f._local_baseline[2]], [f._local_baseline[3], f._local_baseline[4]]
        end
    end
end

@recipe function plot_recipe(curve::Curve{T}, fs::Vector{FWHMmod.FWHM{T}}) where {T<:Number}
    for f in fs
        @series begin
            f
        end
    end
    @series begin
        color --> :blue
        curve
    end
end

@recipe function plot_recipe(crv::Curve{T},
    left::T,
    right::T; subtract_baseline=false,
    local_baseline=false,
    bound=nothing,
    draw_band_centers=false,
    draw_fwhm=false,
) where T
    
    (local_baseline && bound == nothing) && error("You have to provide a bound if local_baseline == true.") |> throw
    (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

    left = T(left)
    right = T(right)
    
    # draw area
    if subtract_baseline
        bh = "end-to-end"
    elseif local_baseline
        bh = "local"
    else
        bh = nothing
    end
    l, r, xl, xr, yl, yr = get_left_right_points(crv.x, crv.y, left, right, bound; baseline_handling=bh)

    if !(local_baseline || subtract_baseline)
        x = [xl; crv.x[l+1:r-1]; xr; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; zero(T); zero(T)]
    elseif subtract_baseline
        x = [xl; crv.x[l+1:r-1]; xr; xl]
        y = [yl; crv.y[l+1:r-1]; yr; yl]
    else
        _, _, ycl, ycr = _endpoint_to_endpoint_baseline(crv.x, crv.y, left, right) # to get extra points located on the curve
        x = [xl; crv.x[l+1:r-1]; xr; xr; xl]
        y = [ycl; crv.y[l+1:r-1]; ycr; yr; yl]
    end

    @series begin
        fillrange := 0
        fillalpha --> 0.5
        fillcolor --> :orange
        linewidth --> 0
        label     --> nothing
        x, y
    end

    if draw_band_centers
        # draw band center
        if subtract_baseline
            bc = band_center(crv, xl, xr, true)
            _, _, yl, yr = _endpoint_to_endpoint_baseline(crv.x, crv.y, xl, xr)
            y0 = lininterp(bc, xl, xr, yl, yr)
        elseif local_baseline
            baseline = baseline_from_points(_local_baseline(crv.x, crv.y, xl, xr, bound)...)
            bc = band_center(crv, xl, xr, baseline)
            _, _, yl, yr = _local_baseline(crv.x, crv.y, xl, xr, bound)
            y0 = lininterp(bc, xl, xr, yl, yr)
        else
            bc = band_center(crv, xl, xr)
            y0 = zero(typeof(crv.y[1]))
        end
        y = lininterp(pmean(bc), crv)
        @series begin
            alpha --> 0.3
            color --> :black
            label --> nothing
            [bc, bc], [y0, y]
        end
    end

    if draw_fwhm
        baseline = local_baseline ? _local_baseline(crv.x, crv.y, xl, xr, bound) : nothing
        @series begin
            fwhm(crv, xl, xr, baseline)
        end
    end
end

function main()
    Random.seed!(42)
    spectrum = NoisySignalIntegration.testdata_1()
    slice_bands = crop(spectrum,  5.0,  40.0)
    slice_noise = crop(spectrum, 40.0, 100.0)

    noise = NoiseSample(slice_noise, 3)
    nm = fit_noise(noise)
    uncertain_spectrum = add_noise(slice_bands, nm)


    # return plot(spectrum, uncertain_spectrum)

    
    position = [15.0, 30.0]
    # widths will fall in the range 2 to 3, with a maximum at 2.5
    width_distribution = scale_shift_beta(2, 2, 3, 4)
    # define a "width bound"
    bds = UncertainBound(position, width_distribution, uncertain_spectrum)
    #@run mc_fwhm(uncertain_spectrum, bds, local_baseline=true)

    # fs = mc_fwhm(uncertain_spectrum, bds; local_baseline=true)
    # return fs
    return uncertain_spectrum, bds

    anim = @animate for k in 1:30
        spectrum_draw_k = get_draw(k, uncertain_spectrum)
        xlims!(12.5, 37.5)
        ylims!(1.2, 2.7)

        xl, xr = get_draw(k, bds[1])
        baseline = _local_baseline(spectrum_draw_k.x, spectrum_draw_k.y, xl, xr, bds[1])
        width1 = fwhm(spectrum_draw_k, xl, xr, baseline)

        xl, xr = get_draw(k, bds[2])
        baseline = _local_baseline(spectrum_draw_k.x, spectrum_draw_k.y, xl, xr, bds[2])
        width2 = fwhm(spectrum_draw_k, xl, xr, baseline)

        plot(spectrum_draw_k, [width1, width2])
    end

    gif(anim, "/tmp/tmp2.gif"; fps=3);
end
