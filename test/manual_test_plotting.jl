## Run this file manually to check that plottings works as expected

using Revise

##
using NoisySignalIntegration


## Imports
using Distributions: MvNormal
using Plots
using Random
using MonteCarloMeasurements

const SAMPLES = 100_000

function get_test_spectrum(seed)
    Random.seed!(seed)
    x = collect(0:0.1:200)
    y = @. exp(-(x-15)^2) + 2 * exp(-(x-30)^2)
    @. y += 1.0 + x*1.5e-2 - (x-50)^3*3e-7
    n = length(x)
    δx = x[2] - x[1] 
    return Curve(x, y .+ (get_cov(δx, n, 0.1, 0.5) |> MvNormal |> rand))
end

spec = get_test_spectrum(1);
noise = NoiseSample(crop(spec, 50, 200), 3)
spec = crop(spec, 10, 40)

nm = fit_noise(noise)
unc_spec = add_noise(spec, nm, SAMPLES)

ubnd1 = UncertainBound(
    Particles(SAMPLES, scale_shift_beta(2.0, 2.0, 10.0, 11.0)),
    Particles(SAMPLES, scale_shift_beta(2.0, 2.0, 19.0, 20.0))
)
ubnd2 = UncertainBound(
    Particles(SAMPLES, scale_shift_beta(2.0, 2.0, 25.0, 26)),
    Particles(SAMPLES, scale_shift_beta(2.0, 2.0, 33.0, 34.0))
)

ubnd3 = UncertainBound(
    15.0,
    scale_shift_beta(2.0, 2.0, 3.0, 4.0),
    unc_spec
)

ubnd4, ubnd5 = UncertainBound(
    [15.0, 30.0],
    scale_shift_beta(2.0, 2.0, 3.0, 4.0),
    unc_spec
)

areas = mc_integrate(unc_spec, [ubnd1, ubnd2, ubnd3, ubnd4, ubnd5])

# test add uncorrelated noise
add_noise(spec, GaussianNoiseModel(1.0)) |> x -> mcplot(x.y, 100, alpha=0.1)

##
plot(unc_spec, [ubnd1, ubnd2]) 
##
plot(unc_spec, ubnd1)

##
plot(spec, unc_spec, 4)

##
plot(noise, nm; draws=2)