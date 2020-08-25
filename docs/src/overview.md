# Package Overview

## Workflow

`NoiseSignalIntegration.jl` estimates the uncertainty in numeric integrals based on the noise level and uncertainty
in placing integration bounds. This is achieved by performing the integration many times while varying the noise and integration bounds.

The package uses a custom datatype [`Curve`](@ref) to represent the xy-data that shall be integrated.
`Curve` wraps two vectors of identical element type and length. It was introduced mainly for convenience
and simpler plotting.

From a `Curve` object, a [`NoiseSample`](@ref) can be derived. A `NoiseSample` is required to determine the noise amplutide
and autocorrelation length (if the noise is strongly correlated, as is often the case in FTIR spectra).
With the noise parameters, a noise model can be constructed. This is either a [`GaussianNoiseModel`](@ref) (uncorrelated Gaussian noise)
or a [`MvGaussianNoiseModel`](@ref) (correlated Gaussian noise).

From the `Curve` object and noise model an [`UncertainCurve`](@ref) object can be constructed. It contains random samples of the original data with varying noise. The `UncertainCurve` object is the first input required for the actual integration function [`mc_integrate`](@ref).

!!! tip "Crop your data to the relevant region"
    While your dataset should contain a somewhat lengthy portion of noise for the noise analysis step,
    you should not include this portion of the data in the actual integration, as this will only
    decrease performance while not offering any benefits. You should always [`crop`](@ref) your data
    to only include the relevant signals you want to integrate (see also [Usage Guide](@ref)).

The second and last input for the integration function is one or several integration bounds.
They are represented by [`UncertainBound`](@ref) objects which are defined from one or two 
probability distributions that encode the uncertainty in the integration bounds.

In the final step, the [`mc_integrate`](@ref) function takes the `UncertainCurve` and `UncertainBound`(s) and integrates
each random sample using the trapezoidal rule. Each `UncertainBound` yields one uncertain area that is returned as a
[Particles](https://baggepinnen.github.io/MonteCarloMeasurements.jl/stable/api/#MonteCarloMeasurements.Particles) object from
the [MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl) package.
Uncertainty propagation with the resulting areas thus works as described in the documentation of `MonteCarloMeasurements.jl`,
merely by performing calculations with `Particles` objects.

!!! tip "Swapping the core integration function"
    You can swap the core integration function that `mc_integrate` uses to something more accurate, if needed.
    To do so, pass your integration function as keyword argument `intfun`.
    Note, that your integration function needs to have the same call signature as [`trapz`](@ref).

    See also: documentation of [`mc_integrate`](@ref).

## Data requirements

Due to the internal workings of the package, the input data needs to fulfill some basic requirements:

!!! warning "Data must be ordered"
    In general, the data that you analyze must be ordered from low to high x-values.
    If your data is not ordered, you should run the `Base.sort()` function on your input data
    (`Curve` object) once in the beginning of your analysis:

    **Example**

    ```jldoctest
    julia> using MCIntegrate

    julia> c = Curve([2, 6, 1], [4, 12, 2]);
    
    julia> c = sort(c)
    Curve{Int64}, 3 datapoints
    (1, 2)
    (2, 4)
    (6, 12)
    ```

!!! warning "x-grid must be uniform when analyzing correlated noise"
    Correlated noise can only be analyzed, if the x-data is evenly spaced.
    If this is not the case for your input data, you can use [Interpolations.jl](http://juliamath.github.io/Interpolations.jl/latest/)
    to interpolate your data on an even x-grid.

    **Example**

    ```@example evengrid
    using MCIntegrate, Interpolations, Plots

    c = testdata_3() # test dataset with uneven grid spacing
    diff(c.x)
    ```

    ```@example evengrid
    δx = minimum(diff(c.x))
    x_even = collect(minimum(c.x):δx:maximum(c.x))
    interp = LinearInterpolation(c.x, c.y)
    c_even = Curve(x_even, interp(x_even))

    plot(c; label="original data")
    plot!(c_even; label="interpolated data")
    ```

- UncertainBound: Width distributions must not have support < 0
