module Temp

using NoisySignalIntegration
using Plots
using Random

function main()
    Random.seed!(42)
    spectrum = NoisySignalIntegration.testdata_1()
    slice_bands = crop(spectrum,  5.0,  40.0)
    slice_noise = crop(spectrum, 40.0, 100.0)

    plot(slice_bands; label="bands")
    plot!(slice_noise; label="noise")

    noise = NoiseSample(slice_noise, 3)
    nm = fit_noise(noise)
    uncertain_spectrum = add_noise(slice_bands, nm)
    position = [15.0, 30.0]
     # widths will fall in the range 2 to 3, with a maximum at 2.5
    width_distribution = scale_shift_beta(2, 2, 3, 4)
    # define a "width bound"
    bds = UncertainBound(position, width_distribution, uncertain_spectrum)
    ctrs = mc_bandcenter(uncertain_spectrum, bds, local_baseline=true)
    @show mc_bandcenter(uncertain_spectrum, bds; local_baseline=true)
    uncertain_spectrum, bds, ctrs
    animate_draws(uncertain_spectrum, bds; filepath="/tmp/tmp.gif", local_baseline=true, draw_band_centers=true)
    plot(uncertain_spectrum, bds; local_baseline=true, draw_band_centers=true)
end

end
