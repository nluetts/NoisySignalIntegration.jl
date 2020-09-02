# Case studies

## Raman spectra

Raman spectra are typically measured with a CCD camera where there is no obviouse autocorrelation of the noise.
The integration works in general as outlined in the [Usage Guide](@ref). The main difference is the analysis 
and generation of noise. Instead of a multivariate `MvGaussianNoiseModel`, we have to use a `GaussianNoiseModel`:

```@example Raman
# simulate a Raman spectrum with two bands

using Plots: plot, plot!, histogram
using Random: seed!
using Statistics
using MonteCarloMeasurements

using NoisySignalIntegration


spectrum = NoisySignalIntegration.testdata_2()
plot(spectrum, label="simulated spectrum")
```

```@example Raman
# crop the bands and the noise sample

bands = crop(spectrum, 10, 40) - 1 # we can directly subtract 1 from the curve to remove the offset
noise = NoiseSample(crop(spectrum, 40, 100), 3)

plot(bands; label="bands")
plot!(noise; label="noise sample")
```

```@example Raman
# build the noise model

nm = noise |> std |> GaussianNoiseModel

uncertain_spectrum = add_noise(bands, nm)

plot(uncertain_spectrum)
```

```@example Raman
# plot some MC draws to check

bounds = UncertainBound([15., 30.], scale_shift_beta(2, 2, 4.5, 5), uncertain_spectrum)

plot(uncertain_spectrum, bounds; subtract_baseline=false)
```

```@example Raman
# integrate the bands and calculate the ratio

area1, area2 = mc_integrate(uncertain_spectrum, bounds; subtract_baseline=false)

histogram(area1/area2)
```

## Mexican hat wavelet

This example is more of a basic sanity check for the integration algorithm.
The integration of a mexican hat wavelet should yield an integral of 0.
For a noisy curve with uncertain integration bounds, the results should still
be compatible with an integral of 0, i.e. the 95% confidence interval of the
resulting area should include the value 0. If we decrease the noise and the
uncertainty in the integration bounds, the result should come closer and closer
to 0.

We go once through the regular workflow. First, we generate the test data:

```@example mh
using NoisySignalIntegration
using Plots
using Random: seed!

seed!(7)

function mexican_hat_curve(noise_model)
    x = -20:0.1:100
    y = @. (1 - x^2)*exp(-x^2/2)
    mh = Curve(x, y)
    uncertain_mh = add_noise(mh, noise_model, 1)
    return NoisySignalIntegration.get_draw(1, uncertain_mh)
end

noisy_mh_curve = mexican_hat_curve(MvGaussianNoiseModel(0.1, 0.05, 0.5))

plot(noisy_mh_curve; label="noisy mexican hat")
```

Analysis of noise:

```@example mh
noise = NoiseSample(crop(noisy_mh_curve, 5, 100))
plot(noise; label="noise sample")
```

```@example mh
noise_model = fit_noise(noise)
```

```@example mh
plotautocovfit(noise, noise_model)
```

Check that generated noise samples look realistic:

```@example mh
plot(noise, noise_model; size=(500, 600))
```

Definition of integration bounds:

```@example mh
uncertain_mh_curve = add_noise(noisy_mh_curve, noise_model)
# definition using distributions for start and end point
bnd = UncertainBound(scale_shift_beta(2, 2, -5, -4), scale_shift_beta(2, 2, 4, 5))
plot(uncertain_mh_curve, bnd; subtract_baseline=false, size=(500, 600), xlim=(-25, 25))
```

Integration:

```@example mh
area = mc_integrate(uncertain_mh_curve, bnd; subtract_baseline=false)

histogram(area; label="area")
```

As expected, the mean of the derived areas is close to zero.

Now we look what happens if we successively decrease the noise level and
uncertainty in the integration window position. We define a function
that covers the complete integration workflow for a certain scaling factor
`f` that controls the noise level and variability of integration windows:

```@example mh
function analyze_noisy_mh(f)
    noisy_mh_curve = mexican_hat_curve(MvGaussianNoiseModel(0.01, 0.05*f, 0.1))
    noise = NoiseSample(crop(noisy_mh_curve, 5, 100))
    noise_model = fit_noise(noise)
    println(noise_model)
    uncertain_mh_curve = add_noise(noisy_mh_curve, noise_model)
    bnd = UncertainBound(scale_shift_beta(2, 2, -4.5-0.5*f, -4.5+0.5*f), scale_shift_beta(2, 2, 4.5-0.5*f, 4.5+0.5*f))
    return mc_integrate(uncertain_mh_curve, bnd; subtract_baseline=false)
end
```

We apply the analysis function to several `f` factors and plot the histograms
of the resulting area distributions:

```@example mh
seed!(2)

areas = [analyze_noisy_mh(f) for f in (1, 0.5, 0.25, 0.05)]

histogram(
    [a.particles for a in areas];
    label=["f = 1" "f = 0.5" "f = 0.25" "f = 0.05"],
    normalize=true,
    alpha=0.3,
    linetype=:stephist,
    fill=true,
    xlim=(-1, 1)
)
```

As we can see, with decreasing uncertainty, the result comes closer and closer to 0.

