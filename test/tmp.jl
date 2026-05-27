using NoisySignalIntegration
using NoisySignalIntegration: Curve, UncertainCurve, UncertainBound
using Random
using Plots

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
    plot(uncertain_spectrum, bds; local_baseline=true, draw_fwhm=true) |> display

    return uncertain_spectrum, bds
end
