# Case studies

## Raman spectra

Raman spectra are typically measured with a CCD camera where there is no obvious autocorrelation of the noise.
The integration works in general as outlined in the [Usage Guide](@ref). The main difference is the analysis 
and generation of noise. Instead of a multivariate `MvGaussianNoiseModel`, we have to use a `GaussianNoiseModel`:

```@setup load_path
push!(LOAD_PATH, "../../src")
```

```@example Raman
# simulate a Raman spectrum with two bands

using Plots: plot, plot!, histogram
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

plot(uncertain_spectrum, bounds; size=(500, 600))
```

```@example Raman
# integrate the bands and calculate the ratio

area1, area2 = mc_integrate(uncertain_spectrum, bounds)

histogram(area1/area2)
```

## Mexican hat wavelet

This example is more of a basic sanity check for the integration algorithm.
The integration of a mexican hat wavelet should yield an integral of 0.
For a noisy curve with uncertain integration bounds, the results should still
be compatible with an integral of 0, i.e. the 95% confidence interval of the
resulting area should include the value 0. If we decrease the noise and the
uncertainty in the integration bounds, the result should come closer and closer
to 0 (the area distribution should become narrower while still covering the 
value 0).

We go once through the regular workflow. First, we generate the test data:

```@example mh
using NoisySignalIntegration
using Plots
using Random: seed!

seed!(7) # seed for reproducibility

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
plot(uncertain_mh_curve, bnd; size=(500, 600), xlim=(-25, 25))
```

Integration:

```@example mh
area = mc_integrate(uncertain_mh_curve, bnd)

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
    return mc_integrate(uncertain_mh_curve, bnd)
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

## [Error propagation](@id example_propagation)

Since integration results are returned as [Particle](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/api/#MonteCarloMeasurements.Particles) objects from the [MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl) package (see [Internals](@ref)), calculating combined uncertainties, for example when deriving abundance ratios from band integrals and calculated intensities, is rather simple. Consider the following spectrum with four bands:

```@example propagation
# simulate an FTIR spectrum with 4 bands

using Plots
using NoisySignalIntegration

spectrum = NoisySignalIntegration.testdata_4()
plot(spectrum, label="simulated spectrum")
```

Let's assume that these bands correspond to 4 different chemical species
(species A, B, C and D from left to right) and we know from quantum chemical
calculations what the intensity of each band should be per unit of substance.
With this information, we can calculate abundance ratios $R_{xy}$ likes so:

```math
R_{xy} = \frac{A_x / I_x}{A_y / I_y}
```

Where $A$ is the band integral and $I$ the calculated intensity for species $x$
or $y$, respectively.

The uncertainty calculation for the abundance ratio is not straightforward if
the distribution of the integrals is asymmetric (which can easily be the case
for larger uncertainties). Using `MonteCarloMeasurement` makes it rather simple,
however.

First, we integrate with `NoisySignalIntegration`:

```@example propagation
# crop spectrum and noise
bands = crop(spectrum, 0, 100)
noise = NoiseSample(crop(spectrum, 100, 200), 3)
# retrieve noise model
nm = fit_noise(noise)
# prepare spectral samples
uspec = add_noise(bands, nm)
# declare integration bounds (several symmetric bands with same width)
bnds = UncertainBound([15., 30., 60., 85.], scale_shift_beta(2, 2, 4, 6), uspec)
# integrate
areas = mc_integrate(uspec, bnds; local_baseline=true)
A_A, A_B, A_C, A_D = areas
```

We inspect the Monte-Carlo draws visually:

```@example propagation
plot(uspec, bnds; size=(500, 600), local_baseline=true)
```

We put in our calculated intensities and their uncertainties (see [`MonteCarloMeasurements.jl` documentation](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/api/#MonteCarloMeasurements.:..-Tuple{Any,Any})):

```@example propagation
using MonteCarloMeasurements

I_A = @samples 10_000   1 .. 1.1  # uniform distribution (.. operator)
I_B = @samples 10_000 1.50 ± 0.2  # normal distribution (± operator)
I_C = @samples 10_000 0.75 ± 0.1
I_D = @samples 10_000   1 .. 1.1
nothing # hide
```

!!! warning "Number of particles must match"
    `MonteCarloMeasurements.jl` propagates uncertainty using a Monte-Carlo
    process with a specific number of samples called "particles". By default,
    `MonteCarloMeasurements.jl` produces 2000 particles when defining an
    uncertain number using the operators `±` and `..` (as of version 0.9.5).
    `NoisySignalIntegration` produces 10000 samples by default. You can use the
    [`@samples`](@ref) macro provided by `NoisySignalIntegration.jl` to increase
    the number of samples produced by the `±` and `..` operators as shown in the
    code example above. The number of samples in all uncertain numbers must
    match when calculating combined uncertainties with
    `MonteCarloMeasurements.jl`, so make sure that this condition is met.

Now we can simply calculate with the retrieved areas and defined intensities and
`MonteCarloMeasurements.jl` takes care of the uncertainty calculation:

```@example propagation
R_BA = (A_B/I_B) / (A_A/I_A)
R_CA = (A_C/I_C) / (A_A/I_A)
R_DA = (A_D/I_D) / (A_A/I_A)
nothing # hide
```

This results in the following, displayed as value ± standard deviation:

```@example propagation
R_BA
```

```@example propagation
R_CA
```

```@example propagation
R_DA
```

We can plot the results as histograms to observe the shape of the distributions:

```@example propagation
plot(
    [R_BA.particles R_CA.particles R_DA.particles],
    label=["B:A" "C:A" "D:A"],
    seriestype=:stephist,
    normalize=true,
    xlabel="abundance ratio",
    ylabel="rel. frequency",
    xlim=(0, 10),
    ylim=(0, 2.5),
    fill=true,
    layout=(3, 1)
)
```

Clearly, the resulting distributions are asymmetric and non-Gaussian, so the
standard deviations do not inform about the level of confidence. You can use
[`StatsBase.jl`](https://juliastats.org/StatsBase.jl/stable/) to calculate
percentiles and confidence intervals:

```@example propagation
using StatsBase: percentile

[percentile(R_BA.particles, p) for p in (2.5, 50, 97.5)]
```

A visual representation like the histogram may better convey the range of uncertainty, though.

Alternatives are a box plot:

```@example propagation
using StatsPlots: boxplot

let
    rep(str) = repeat([str], length(R_BA.particles)) 
    x = [rep("B:A") rep("C:A") rep("D:A")]
    y = [R_BA.particles R_CA.particles R_DA.particles]
    boxplot(x, y, ylabel="abundance ratio", label=nothing)
end
```

Or a violin plot:

```@example propagation
using StatsPlots: violin

let
    rep(str) = repeat([str], length(R_BA.particles)) 
    x = [rep("B:A") rep("C:A") rep("D:A")]
    y = [R_BA.particles R_CA.particles R_DA.particles]
    violin(x, y, ylabel="abundance ratio", label=nothing)
end
```
