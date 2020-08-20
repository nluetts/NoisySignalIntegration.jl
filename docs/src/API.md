# API Reference

```@contents
Pages = ["API.md"]
Depth = 3
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

add_noise
```

### Noise analysis

```@docs
fit_noise

plot_autocov

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