
# Usage guide

As an usage example, we will go through the analysis of FTIR and Raman band signals.

Suppose we measured an FTIR spectrum that looks like the following simulation:

```@example FTIR
using Distributions: MvNormal
using Plots: plot, plot!
using Random: seed!

using MCIntegrate

seed!(42)

x = collect(0:0.1:100)

# simulate FTIR spectrum with two bands
# the true intensity ratio is 1 : 2
bands = @. exp(-(x-15)^2) + 2 * exp(-(x-30)^2)
baseline = @. 1.0 + x*1.5e-2 - (x-50)^3*3e-6
noise = begin
    n = length(x)
    δx = x[2] - x[1] 
    get_cov(δx, n, 0.1, 0.5) |> MvNormal |> rand
end

spectrum = Curve(x, baseline + bands + noise)

plot(spectrum, label="simulated spectrum")
```

In order two apply the MCIntegrate uncertainty analysis, we must perform 4 basic steps:
1. Crop from the spectrum the region that contains the signals and the region that contains a representative sample of the noise
1. Characterize the noise (to be able to simulate it in the Monte-Carlo draws) 
1. Set integration bounds and their associated uncertainties
1. Run the `mcintegrate()` function

## Cropping the spectrum

Let's start by dividing the spectrum into the bands we want to integrate and the noise we
want to analyse. We can do this by using the `crop()` function.

```@example FTIR
slice_bands = crop(spectrum,  5.0,  40.0)
slice_noise = crop(spectrum, 40.0, 100.0)

plot(slice_bands; label="bands")
plot!(slice_noise; label="noise")
```

## Noise analysis

The spectrum has a quite considerable baseline which constitutes a problem when analysing the noise. To prepare the noise spectrum `slice_noise` for analysis, we create a `NoiseSample` object. Upon construction of the `NoiseSample` object, a polynomial is fitted and subtracted to remove the baseline:

```@example FTIR
noise = NoiseSample(slice_noise, 3) # 3 = remove third order polynomial baseline
nothing # hide
```

Plotting the original slice and the `NoiseSample` object shows the data after baseline removal:

```@example FTIR
plot(slice_noise, label="cropped noise")
plot!(noise, label="NoiseSample object")
```

In order to simulate the noise, we must determine its characteristics.
A model is retrieved by fitting the estimated autocovariance:

```@example FTIR
nm = fit_noise(noise)
# plot the fitting result:
plot_autocov(noise, nm)
```

Plotting the model next to the noise object is an important sanity check to verify that the fitting yielded a sensible estimate and that generated noise samples do mimic the experimental noise. They keyword `draws` controls how many random generated noise draws are plotted.

```@example FTIR
plot(noise, nm)
```

## Preparing the spectrum for integration

Now that we have a noise model, we can generate an `UncertainCurve`. An `UncertainCurve` holds random draws of the original spectrum plus noise:

```@example FTIR
uncertain_spectrum = add_noise(slice_bands, nm, 50_000)
nothing # hide
```

If we plot the `uncertain_spectrum`, we get a ribbon plot showing a 95% confidence band:

```@example FTIR
plot(uncertain_spectrum)
```

We can plot also single draws by using the `mcplot()` function from `MonteCarloMeasurements.jl`:

```@example FTIR
using MonteCarloMeasurements

mcplot(uncertain_spectrum; draws=20)
```


## Integration bounds

MCIntegrate deals with uncertainty in placing integration bounds by expressing each bound by one ore more probability distributions. Any continuous, univariate distribution from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) can be used to define integration bounds.

To define integration bounds, the `UncertainBound` object is used.
There are several options available to create an `UncertainBound`:

1. Passing two distributions
1. Passing a position, a distribution, and an `UncertainCurve`
1. Passing a vector of positions, a distribution, and an `UncertainCurve`

### Defining the bound using a start and end point

If two distributions are passed to `UncertainBound()`, they will be used to 

A `LeftRightBound` is used to specify an integration window by a start and an end point. Start and end point are expressed by probability distributions.

Here is an example. Let's say we want to integrate the right peak of our simulated spectrum:

```@example FTIR
plot(crop(spectrum, 20, 40), label="right peak")
plot!([27, 32], [1.3, 1.3]; markershape=:cross, label="integration interval")
```

It looks like integrating from about 27 to 32 would be appropriate, but there is some doubt of the exact location of the integration bounds. Perhaps a reasonable estimate is that the left bound falls in the range from 26 to 27 and the right bound in the range from 32 to 33. This would be expressed with a `LeftRightBound` that is defined using two uniform distributions:

```@example FTIR
using Distributions: Uniform

lrb = LeftRightBound(Uniform(26, 27), Uniform(32, 33))
nothing # hide
```

We can draw random samples from the `LeftRightBound` object and plot a histogram of the samples to visualize the generated integration bounds.

```@example FTIR
using Plots: histogram, histogram!

histogram(sample(lrb, 10_000); label=false)
```

The uniform distribution is of course a bit of an awkward choice, because its probability density suddenly drops to 0, which perhaps does not model one's belief about the position of the integration bound very well. Ont the other hand, the normal distribution is often a natural choice when dealing with uncertainties:

```@example FTIR
using Distributions: Normal

lrb_normal = LeftRightBound(Normal(26.5, 0.5), Normal(32.5, 0.5))

histogram(sample(lrb_normal, 10_000); label=false)
```

However, in this particular case of describing the uncertainty of integration bounds, the tails of the normal distribution are problematic, because they lead to occasional extreme values of the intergration interval, which would not seem realistic.

A compromise between the uniform and normal distribution, that can be used to model the unertainty in the bound start and end values, is a scaled and shifted beta(2, 2) distribution. Its shape resembles the shape of the normal distribution but it is missing the tails. Since a scaled and shifted beta distribution does not ship with Distributions.jl, MCIntegrate includes the function `scaled_shifted_beta(α, β, a, b) ` which can be used to generate a beta(α, β) distribution that has a support region in the interval a to b.

Again, a demonstration may help to explain. We keep the normal distribution for the right bound so we can compare the distributions easily:

```@example FTIR
using Distributions: Normal

lrb_beta_normal = LeftRightBound(scale_shift_beta(2, 2, 26, 27), Normal(32.5, 0.5))

histogram(sample(lrb_beta_normal, 10_000); label=false, normalize=true)
```

### `WidthBound` and `WidthBoundClone`

The `WidthBound` describes the integration window by a location (a single number) and a width (a distribution). This is not merely for convenience. The integration window of a `WidthBound` will always be symmetric around the peak maximum that falls into the range of the location parameter ± the mean of the distribution associated with the width. This is particularly useful if the noise level is rather high and the peak position may move considerable in each Monte-Carlo iteration.
The `WidthBound` "follows" the peak, so to speak.

A `WidthBoundClone` can be generated from a parent `WidthBound`. Clones can be used to integrate further bands with the same window width as the parent in each Monte-Carlo iteration. This is useful if one can assume from a physical argument that the peaks should have the same width.

For example, to define such bounds for the peaks in the simulated spectrum, we would do:

```@example FTIR

wb_1 = WidthBound(15, Uniform(2, 3))
wb_2 = WidthBoundClone(30, wb_1)

histogram( sample(wb_1, slice_bands, 100); label=false, normalize=true)
histogram!(sample(wb_2, slice_bands); label=false, normalize=true)
```

Note that sampling from a `WidthBound` (or its clone) yields samples of integration start and end points, which is the same behavior as in the case of sampling from a `LeftRightBound`. The spectrum has to be passed to the `sample` function to determine the peak position around which the symmetric integration windows shall be placed.

In the above example, the sample size was deliberately choosen rather small, so one can see that the samples of the clone exactly follow the samples of the parent.

## Plotting Monte-Carlo draws

To verify that the integration windows and derived integrals are sensible, it is a good idea to plot a few draws and integrals before running the full Monte-Carlo algorithm. We can do so by passing a spectrum, an array of bounds, and a noise model to the plot function:

```@example FTIR
plot(slice_bands, [wb_1, wb_2], nm; samples=4)
```

We can see from the plot that our estimate for the width of the peaks was perhaps a bit too small, so we retry:

```@example FTIR
wb_1 = WidthBound(15, Uniform(3, 4))
wb_2 = WidthBoundClone(30, wb_1)

plot(slice_bands, [wb_1, wb_2], nm; samples=4)
```

## Running the integration algorithm

The integration is performed with the function `mcintegrate`. We have to pass in the spectrum, integration bounds, noise model, and number of draws:

```@example FTIR
integral_samples = mc_integrate(slice_bands, [wb_1, wb_2], nm; N=10_000)
```

Now we can look at the histogram of the integrals, or of the peak area ration:

```@example FTIR
histogram(integral_samples)
```

```@example FTIR
ratio_samples = integral_samples[:, 2] ./ integral_samples[:, 1]
histogram(ratio_samples)
```

We see that the histogram of peak area ratio peaks around 2, which is what we put into the simulation of the spectrum. Considering the noise and the uncertainty in the integration bounds, we end up with a uncertainty interval of roughly 1 to 3

```@example FTIR
using StatsBase: mean, std, percentile

mean(ratio_samples) |> println
std(ratio_samples) |> println
percentile(ratio_samples, 2.5) |> println
percentile(ratio_samples, 50) |> println
percentile(ratio_samples, 97.5) |> println
```