# NoisySignalIntegration.jl

*NoisySignalIntegration.jl -- A tool to determine uncertainty in numeric
integrals of noisy x-y data.*

`NoisySignalIntegration` implements a method to determine the uncertainty in
numeric integrals of noisy x-y data on the basis of a Monte-Carlo process.  It
can include uncertainty due to noise, baseline subtraction, and placement in
integration bounds.  To do this, the integration is repeated many times while
the noise of the data, baseline, and integration bounds are varied based on a
noise model and user supplied probability distributions.

A predecessor of this package was originally intended to estimate uncertainties
of band signals in FTIR spectra (see [G. Karir et al.,
2019](https://doi.org/10.1039/C9CP00435A)), which is reflected in the example
given in the [Usage Guide](@ref).

**Table of Contents**

```@contents
Pages = ["overview.md", "guide.md", "examples.md", "baseline.md", "API.md", "internals.md"]
Depth = 4
```
