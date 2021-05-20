# Usage Guide

As a more detailed usage example, we will go through the analysis of a simulated
FTIR spectrum.

!!! info "A note on plotting"
    `NoisySignalIntegration` provides several
    "[recipes](http://docs.juliaplots.org/latest/recipes/)" for
    [`Plots.jl`](http://docs.juliaplots.org/latest/) to easily plot the various
    (interim) results. Often, merely calling `plot()` and passing in data types
    from `NoisySignalIntegration` will work. Examples are included in this
    guide.

Suppose our spectrum looks like the following simulation:

```@example FTIR
using Distributions, NoisySignalIntegration, Plots
using Random: seed!

spectrum = NoisySignalIntegration.testdata_1()
plot(spectrum, label="simulated spectrum")
```

In order two apply the `NoisySignalIntegration` uncertainty analysis, we must perform 4 basic steps:

1. From the spectrum, crop the region that contains the signals and the region
   that contains a representative sample of the noise
1. Characterize the noise (to be able to simulate it in the Monte-Carlo draws) 
1. Set integration bounds and their associated uncertainties
1. Run the [`mc_integrate`](@ref) function

## [Cropping the spectrum](@id crop)

Let's start by dividing the spectrum into the bands we want to integrate and the
noise we want to analyse. We can do this by using the [`crop`](@ref) function.

```@example FTIR
slice_bands = crop(spectrum,  5.0,  40.0)
slice_noise = crop(spectrum, 40.0, 100.0)

plot(slice_bands; label="bands")
plot!(slice_noise; label="noise")
```

## Noise analysis

The spectrum has a quite considerable baseline which constitutes a problem when
analysing the noise. To prepare the noise spectrum `slice_noise` for analysis,
we create a [`NoiseSample`](@ref) object. Upon construction of the `NoiseSample`
object, a polynomial is fitted and subtracted to remove the baseline:

```@example FTIR
noise = NoiseSample(slice_noise, 3) # 3 = remove third order polynomial baseline
nothing # hide
```

Plotting the original slice and the `NoiseSample` object shows the data after
baseline removal:

```@example FTIR
plot(slice_noise, label="cropped noise")
plot!(noise, label="NoiseSample object")
```

In order to simulate the noise, we must determine its characteristics. A model
is retrieved by fitting the estimated autocovariance:

```@example FTIR
nm = fit_noise(noise)

# plot the fitting result:
plotautocovfit(noise, nm);
# create zoomed inset:
lens!([0, 1.5], [-1e-3, 3e-3], inset = (1, bbox(0.3, 0.3, 0.3, 0.3)))
```

Plotting the model next to the noise object is an important sanity check to
verify that the fitting yielded a sensible estimate and that generated noise
samples do mimic the experimental noise.

```@example FTIR
plot(noise, nm)
```

!!! info "plotting more samples"
    They keyword `draws` controls how many random generated noise draws are
    plotted:
    
    ```julia
    plot(noise, nm; draws=5) # draw 5 instead of 3 (default) noise draws
    ```

    This also works when you [plot Monte-Carlo draws](@ref plot_mc_draws).

## Preparing the spectrum for integration

Now that we have a noise model, we can generate an [`UncertainCurve`](@ref).
An `UncertainCurve` holds random draws of the original spectrum plus noise:

```@example FTIR
uncertain_spectrum = add_noise(slice_bands, nm, 50_000) # generate 50_000 random samples
nothing # hide
```

If we plot the `uncertain_spectrum`, we get a ribbon plot showing a 95%
confidence band:

```@example FTIR
plot(uncertain_spectrum)
```

We can also plot single draws by using the `mcplot()` function from
`MonteCarloMeasurements.jl`:

```@example FTIR
using MonteCarloMeasurements

mcplot(uncertain_spectrum; draws=20)
```


## [Integration bounds](@id bounds_guide)

`NoisySignalIntegration` deals with uncertainty in placing integration bounds by
expressing each bound by one ore more probability distributions. Any continuous,
univariate distribution from
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl) can be used
to define integration bounds.

To define integration bounds, the `UncertainBound` object is used. There are
several options available to create an `UncertainBound`:

1. Passing two distributions
1. Passing a position, a distribution, and an `UncertainCurve`
1. Passing a vector of positions, a distribution, and an `UncertainCurve`

### Defining an `UncertainBound` using a start and end point

If two distributions are passed to `UncertainBound()`, they will be interpreted
as start and end points for the integration window with the uncertainty of these
points being expressed by the spread of the distributions.

For example, let's say we want to integrate the right peak of our simulated
spectrum:

```@example FTIR
plot(crop(spectrum, 20, 40), label="right peak")
plot!([27, 32], [1.3, 1.3]; markershape=:cross, label="integration interval")
```

It looks like integrating from about 27 to 32 would be appropriate, but there is
some doubt of the exact location of the integration bounds. Perhaps a reasonable
estimate is that the left bound falls in the range from 26 to 27 and the right
bound in the range from 32 to 33. This would be expressed with a
`UncertainBound` that is defined using two uniform distributions:

```@example FTIR
lrb = UncertainBound(Uniform(26, 27), Uniform(32, 33)) # 10 k samples by default
nothing # hide
```

Upon creation of the `UncertainBound` object, pseudo random samples of the
integration start and end point are drawn. If we do not provide the number of
samples, it will default to 10 000. We can plot the bound as a histogram to see
the distribution of the start and end point:

```@example FTIR
histogram(lrb; label=["start" "end"])
```

The uniform distribution is of course a bit of an awkward choice, because its
probability density suddenly drops to 0, which perhaps does not model one's
belief about the position of the integration bounds very well.

Due to the central limit theorem and the general applicability of the normal
distribution, it is often a natural choice when dealing with uncertainties:

```@example FTIR
lrb_normal = UncertainBound(Normal(26.5, 0.5), Normal(32.5, 0.5), 12_000) # we draw 12_000 samples, just to illustrate how it works

histogram(lrb_normal; label=["start" "end"])
```

However, in this particular case of describing the uncertainty of integration
bounds, the tails of the normal distribution are problematic, because they lead
to occasional extreme values of the integration interval, which would not seem
realistic.

A compromise between the uniform and normal distribution is a scaled and shifted
beta(2, 2) distribution. Its shape resembles the shape of the normal
distribution but it is missing the tails. Since a scaled and shifted beta
distribution does not ship with
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl),
`NoisySignalIntegration` includes the function [`scale_shift_beta`](@ref)`(α, β,
a, b)` which can be used to generate a beta(α, β) distribution that has a
support in the interval `a` to `b`.

Again, a demonstration may help to explain (we keep the normal distribution for
the right bound so we can compare the distributions):

```@example FTIR
lrb_beta_normal = UncertainBound(scale_shift_beta(2, 2, 26, 27), Normal(32.5, 0.5))

histogram(lrb_beta_normal; label=["start" "end"], normalize=true)
```

### Defining an `UncertainBound` using a position, width, and `UncertainCurve` (symmetric bands)

For some spectra, one can assume that a band is more or less symmetric. In this
case, it may be better to define an integration window not by a start and end
point but merely by a width, and integrate the band symmetrically around its
peak position ± half this width.

To accomplish this, one has to construct an `UncertainBound` object by passing a
`position` (the peak position of the symmetic band), a distribution that
describes the uncertainty in the integration window width (here
`width_distribution`), and an `UncertainCurve` that holds samples of the
spectrum (here `uncertain_spectrum`):

```@example FTIR
position = 15.0
 # widths will fall in the range 2 to 3, with a maximum at 2.5
width_distribution = scale_shift_beta(2, 2, 2, 3)
# define a "width bound"
wb = UncertainBound(position, width_distribution, uncertain_spectrum)
nothing # hide
```

From the provided data, the `UncertainBound` object is created as follows:
- The width distribution is used to draw as many random samples of the
  integration window width $w$ as the `uncertain_spectrum` contains spectral
  samples
- For each spectral sample in `uncertain_spectrum`, the peak position $p_x$ in
  the range `position` ± $\frac{w}{2}$ is retrieved
- The peak position $p_x$ is used to define the start and end point of the
  integration window for each spectral sample, $p_x - \frac{w}{2}$ and $p_x +
  \frac{w}{2}$

Therefore, what is stored in the `UncertainBound` object are again start and end
points for the integration. We can verify this by plotting another histogram:

```@example FTIR
histogram(wb; label=["start" "end"], normalize=true, linewidth=0)
```

The crucial difference compared to the bound defined from two distributions is
that the start and end points are now placed symmetrically with respect to the
band's peak position. The `UncertainBound` now "follows" the peak in each
Monte-Carlo draw, so to speak.

### Defining an `UncertainBound` using several positions, a width, and `UncertainCurve` (several symmetric bands with same width)

If, for example from a physical argument, we can say that two bands should have
the same width, we can constrain our `UncertainBound`s even further: we can
create several bounds that share the exact same integration window width in each
draw.

All we have to do is to provide the constructor of `UncertainBound` not with a
single position, but with an array of several positions:

```@example FTIR
positions = [15.0, 30.0]

wb_1, wb_2 = UncertainBound(positions, width_distribution, uncertain_spectrum)
nothing # hide
```

Note that the constructor will then return an array of `UncertainBound` objects
which we unpacked into the variables `wb_1` and `wb_2` in the example above.

The histograms of the start and end points looks like this:

```@example FTIR
histogram( wb_1; label=["start 1" "end 1"], normalize=true, linewidth=0)
histogram!(wb_2; label=["start 2" "end 2"], normalize=true, linewidth=0)
```

It is not obvious from the histograms that the widths of the integration windows
stored in `wb_1` and `wb_2` are identical, so we calculate and print them here
manually to prove this:

```@example FTIR
for i in 1:10
    l1 = wb_1.left.particles[i]
    l2 = wb_2.left.particles[i]
    r1 = wb_1.right.particles[i]
    r2 = wb_2.right.particles[i]

    println("draw $i: width 1 = ", r1 - l1, " width 2 = ", r2 - l2)
end
```

!!! warning "Watch out for the support of your width distribution"
    Note that the distribution that you pass to `UncertainBound` along with a
    position/positions must not allow for negative values (i.e. its support must
    end before 0). Keep in mind that a normal distribution, for
    example, has support from -∞ to ∞, so it is a poor choice here.

## [Plotting Monte-Carlo draws](@id plot_mc_draws)

To verify that the integration windows and derived integrals are sensible, it is
a good idea to plot a few draws and integrals before running the full
Monte-Carlo algorithm. We can do so by passing an `UncertainCurve` and an array
of `UncertainBound`s to the plot function:

```@example FTIR
plot(uncertain_spectrum, [wb_1, wb_2]; size=(400, 500), local_baseline=true)
```

We can see from the plot that our estimate for the width of the peaks was
perhaps a bit too small, so we retry:

```@example FTIR
width_distribution = scale_shift_beta(2, 2, 3, 4) # width will fall in the range [3, 4]
wb_1, wb_2 = UncertainBound(positions, width_distribution, uncertain_spectrum)

plot(uncertain_spectrum, [wb_1, wb_2]; size=(400, 500), local_baseline=true)
```

## Running the integration algorithm

The integration is performed with the function [`mc_integrate`](@ref). We have
to pass in the uncertain spectrum and integration bounds. Since we pass in two
integration bounds, we retrieve two areas:

```@example FTIR
area_1, area_2 = mc_integrate(uncertain_spectrum, [wb_1, wb_2]; local_baseline=true)
```

We can look at the histogram of the integrals:

```@example FTIR
histogram([area_1.particles, area_2.particles]; label=["band area 1" "band area 2"])
```

Or of the peak area ratio, simply by calculating with the retrieved areas:

```@example FTIR
ratio = area_1 / area_2
histogram(ratio.particles; label="band area ratio (band 1/band 2)")
```

We see that the histogram of peak area ratio peaks around 0.5, which is what we
put into the simulation of the spectrum.

We can use some basic statistical functions to characterize the result:

```@example FTIR
using StatsBase: mean, std, percentile

mean(ratio)
```

```@example FTIR
std(ratio)
```

```@example FTIR
percentile(ratio, 2.5)
```

```@example FTIR
percentile(ratio, 50)
```

```@example FTIR
percentile(ratio, 97.5)
```

We find that, considering the noise and the uncertainty in the integration
bounds, we end up with a 95% uncertainty interval of roughly 0.4 to 0.7.

!!! info "Sensitivity analysis"
    It is perfectly valid to create several `UncertainBound`s for *one and the
    same band* and feed them into `mc_integrate()`, e.g. to perform a sensitivity
    analysis on how much the result depends on the kind and parameters of the
    bounds.

You find more usage examples on the next page. In particular, check out the
[error propagation example](@ref example_propagation) to see how to proceed with
uncertainty calculations with the retrieved areas using
[`MonteCarloMeasurements.jl`](https://github.com/baggepinnen/MonteCarloMeasurements.jl).