using NoisySignalIntegration
using NoisySignalIntegration: get_draw, _local_baseline, lininterp
using MonteCarloMeasurements
using Plots
using RecipesBase
using Random
using Debugger

module FWHMmod

struct FWHM{T}
    full_width::T
    _local_baseline::Union{Nothing, Tuple{T, T, T, T}}
    _peak_position::T
    _start_position::T
    _half_maximum::T
    _half_maximum_offset::T
end

struct UncertainFWHM{T}
    samples::Vector{FWHM{T}}
end

function value(uf::UncertainFWHM{T}) :: MonteCarloMeasurements.Particles{T, N} where {T,N}
    particles([f.full_width for f in uf.samples])
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
    return [FWHMmod.UncertainFWHM(ws) for ws in eachcol(widths)]
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
        # label = "$(([yi, hm, yj] .|> x -> round(x, digits=3)))"
        # println(label)
        # plot(xs, ys, ylims=(1, 2.1), label=label)
        # plot!([left, right], [hm + ybx, hm + ybx])
        # scatter!([xs[i], xs[i + 1]], [yi + ybx, yj + ybx])
        # frame(anim)
        if yj < hm <= yi 
            # find x-value that corresponds best to max/2
            x_fwhm_right = lininterp(hm, yi, yj, xs[i], xs[i + 1])
            # scatter!([x_fwhm_right], [hm + ybx])
            # frame(anim)
            break
        end
    end

    # gif(anim, "/tmp/tmp.gif", fps=25, loop=1)

    # p = plot(xs, ys)
    # plot!(p, [x_fwhm_left, x_fwhm_right], [hm + ybx, hm + ybx])
    # display(p)

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

@recipe function plot_recipe(
    uc::UncertainCurve{T, N},
    fwhms::Vector{FWHMmod.FWHM{T}}
    ;
    draws=3,
    subtract_baseline=true
) where {T, N}

    legend := :none
    layout := (draws + 1, 1)
    link := :both
    size --> (500, 600)
    
    mean_uc = mean(uc)
    
    for i ∈ 0:draws
        for (j, b) in enumerate(bnds)
            @series begin
                fillcolor := j % 2 == 1 ? :red : :orange
                subplot := i + 1
                bound := b
                if i == 0
                    mean_uc, mean(b)...
                else
                    get_draw(i, uc), get_draw(i, b)...
                end
            end
        end
        @series begin
            subplot := i + 1
            if i == 0
                # mean spectrum
                seriescolor := SECONDARY_COLOR
                yguide := "mean"
                mean_uc
            else
                seriescolor := PRIMARY_COLOR
                yguide := "sample $(i)"
                get_draw(i, uc)
            end
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
