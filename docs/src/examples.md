# Case studies

## Raman spectra

Raman spectra are in general measured with a CCD camera where there is no obviouse autocorrelation of the noise.
The integration works in general as outlined in the [Usage Guide](@ref). The main difference is the analysis 
and generation of noise. Instead of a multivariate `MvGaussianNoiseModel`, we have to use a `GaussianNoiseModel`:


```@example 
using Plots: plot, plot!, histogram
using Random: seed!
using Statistics
using MonteCarloMeasurements

using MCIntegrate

seed!(42)

function get_spectrum(σ=0.1)
    x = collect(0:0.1:100)
    # simulate Raman spectrum with two bands
    # the true intensity ratio is 1 : 2
    bands = @. exp(-(x-15)^2) + 2 * exp(-(x-30)^2)
    baseline = @. 1.0 + 1.5e-4x^2 - 3e-6x^3
    noise = randn(length(x)) * σ
    Curve(x, baseline + bands + noise)
end

spectrum = get_spectrum(0.001)

plot(spectrum, label="simulated spectrum")

bands = crop(spectrum, 10, 40) - 1 # remove offset
noise = NoiseSample(crop(spectrum, 40, 100), 3)

plot(bands; label="bands")
plot!(noise; label="noise sample")

std(noise)

nm = noise |> std |> GaussianNoiseModel

us = add_noise(bands, nm)

plot(us)

mcplot(us)

bnds = UncertainBound([15., 30.], scale_shift_beta(2, 2, 4.5, 5.5), us)

plot(us, bnds)

area1, area2 = mc_integrate(us, bnds)

histogram(area1/area2)
```