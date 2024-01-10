# Baseline Handling

There are several ways to handle baseline correction when working with the package.
The easiest method is to use the build-in local baseline correction which
assumes a linear baseline between the start and end point of the integration
window. It is envoked by using the keyword argument `subtract_baseline`
(deprecated in v0.2) or `local_baseline`. The difference of these methods is discussed
below.

Otherwise, one can subtract a baseline from the data in a preprocessing step.
This can be done either before or after generating an [`UncertainCurve`](@ref).
If a baseline is subtracted from the `UncertainCurve`, it is possible to
account for uncertainty in the baseline correction, e.g. by subtracting
baselines generated using a Gaussian process.

## Build-in

Local linear baseline subtraction can be achieved by passing
`local_baseline=true` to the [`mc_integrate`](@ref) function. To visualize the
integrated area, the same keyword argument can be passed to the `plot()`
function when [plotting Monte-Carlo draws](@ref plot_mc_draws).

The keyword argument `subtract_baseline=true` is also supported, but its use is
deprecated. The difference of `local_baseline` and `subtract_baseline` can be
visualized when animating draws of curves and integration bound samples (using
`Plots.@animate`): 

```@setup load_path
push!(LOAD_PATH, "../../src")
```

```@eval
using NoisySignalIntegration
using Plots
using Random: seed!

seed!(1)
nsi = NoisySignalIntegration

let
    n = 20
    c = crop(nsi.testdata_1(), 0, 50)
    uc = add_noise(c, GaussianNoiseModel(0.1))
    ubleft = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 4.0, 5.0), uc)
    ubright = UncertainBound(30.0, scale_shift_beta(2.0, 2.0, 6.0, 7.0), uc)

    spany = [
        mm(curve.y for curve ∈ [nsi.get_draw(i, uc) for i ∈ 1:n]) |> mm
        for mm in (minimum, maximum)
    ]
    spany = (spany[1]*0.9, spany[2]*1.1)
    anim = @animate for i in 1:n
        kw = Dict(:ylim => spany, :legend => :topleft)
        p1 = plot(nsi.get_draw(i, uc), [ubleft, ubright], i; subtract_baseline=true, label=["subtract_baseline" "" ""], kw...)
        p2 = plot(nsi.get_draw(i, uc), [ubleft, ubright], i; local_baseline=true, label=["local_baseline" "" ""], kw...)
        plot(p1, p2; layout=(2, 1), size=(400, 300))
    end
    gif(anim, "baseline_anim.gif", fps=5)
    nothing
end
```

![build-in baseline handling animation](baseline_anim.gif)

As the figure above shows, the baseline varies considerably more when
`subtract_baseline` is used, compared to `local_baseline`. The former uses the
exact start and end points of the integration window in a particular draw while
the latter uses the start and end point distributions to determine a weighted
average for y-values at integration bounds of a particular draw. Especially at high noise
levels, rather extreme local baseline estimates can follow from using
`subtract_baseline` and the overall uncertainty may be overestimated. It is
thus recommended to use the `local_baseline` method.

## Preprocessing

### Simple baseline subtraction

The data can simply be preprocessed before a [`Curve`](@ref) object is created.
For example, one can mask the signals that shall be integrated (using a filter
or calls to the [`crop`](@ref) and [`stitch`](@ref) functions) and fit a
[polynomial](https://github.com/JuliaMath/Polynomials.jl) or
[smoothing spline](https://github.com/nignatiadis/SmoothingSplines.jl) to be
subtracted:



```@setup simple_baseline
using NoisySignalIntegration, Plots
crv = NoisySignalIntegration.testdata_1()
```

```@example simple_baseline
using Polynomials, SmoothingSplines

polyfit, splinefit = let
    no_signals = stitch(crop(crv, 0, 12), crop(crv, 32, 300)) # remove signals
    x = no_signals.x
    y = no_signals.y
    pfit = Polynomials.fit(x, y, 5)
    sfit = SmoothingSplines.fit(SmoothingSpline, x, y, 2000.0)
    pfit, sfit
end

crv_baseline_corrected = crv - predict(splinefit, crv.x)

plot(crv; label="data", alpha=0.3, color=:black, legend=:outertopright, size=(800, 400))
plot!(polyfit, 0, 100; label="polynomial baseline fit")
plot!(crv.x, predict(splinefit, crv.x); label="spline baseline fit")
plot!(crv_baseline_corrected; label="data - spline baseline")

```

### Uncertain baseline

It may be desirable to account for uncertainty in the subtracted baseline.
One way to achieve this is to fit a [Gaussian process](https://github.com/STOR-i/GaussianProcesses.jl) to the dataset while
excluding signals. For example, consider a dataset with a relatively broad
band:

```@setup ubaseline
using NoisySignalIntegration, Plots, Polynomials

baseline = Polynomial([1, -0.003, 1e-5, 3e-7]).(-150:0.5:150)
nm = MvGaussianNoiseModel(0.5, 0.1, 1.5);

crv = NoisySignalIntegration.generate_testdata(
    0:0.5:300,
    [(15, 40, 5), (30, 140, 10), (10, 150, 15), (5, 155, 8), (25, 170, 10)],
    nm;
    seedvalue=42,
    baseline=baseline
)
```

```@example ubaseline
plot(crv, label=nothing)
```

A Gaussian process can be used to approximate the baseline, while allowing for
higher uncertainty in the regions where the course of the baseline is masked by
signals:


```@example ubaseline
using GaussianProcesses, MonteCarloMeasurements

gp = let
    # cut out signals
    no_signals = stitch(crop(crv, 0, 25), crop(crv, 55, 110), crop(crv, 200, 300))
    # fit Gaussian process
    GP(no_signals.x, no_signals.y, MeanZero(), SE(4.0, 0.1))
end

# draw random baselines ...
fs = rand(gp, crv.x, 250)

# ... and plot
plot(crv.x, fs; alpha=0.05, color=:red, dashed=:dashed, legend=nothing)
plot!(crv)
```

The "amplitude" of the baseline's uncertainty in the range of the signals
can be tweaked by the hyperparameters of the Gaussian process, in particular
the autocorrelation length scale of the covariance function (in the example
`SE()`, the "squared exponential" function, i.e. Gaussian covariance function,
with a length scale of 4.0 units).


The baseline samples can be subtracted from an `UncertainCurve` (here `ucrv`) as follows:

```julia
# make sure the number of samples matches, here 100 000 samples
ubaseline = rand(gp, crv.x, 100_000) |> transpose |> collect |> Particles
ucrv_baseline_corrected = UncertainCurve(ucrv.x, ucrv.y - ubaseline)
```

The resulting `UncertainCurve` now includes not only uncertainty due to noise
but also due to the baseline correction.
