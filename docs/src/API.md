# API Reference

```@contents
Pages = ["API.md"]
Depth = 4
```

## Types

```@docs
Curve

NoiseSample

UncertainCurve

GaussianNoiseModel

MvGaussianNoiseModel

UncertainBound
```

## Functions

### Manipulation of Curves

```@docs
crop

stitch

add_noise
```

### Noise analysis

```@docs
fit_noise

plotautocovfit

```

### Statistics

```@docs
scale_shift_beta
```

### Integration

```@docs
mc_integrate

trapz
```

### Band Centers

```@docs
mc_bandcenter
````

## Macros

```@docs
@samples
```

## Misc

```@docs
animate_draws
```
