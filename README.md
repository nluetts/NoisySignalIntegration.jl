# NoisySignalIntegration.jl

*A tool to determine uncertainty in numeric integrals of noisy x-y data.*

`NoisySignalIntegration` implements a method to determine the uncertainty in
numeric integrals of noisy x-y data on the basis of a Monte-Carlo process.  It
can include uncertainty due to noise, baseline subtraction, and placement in
integration bounds.  To do this, the integration is repeated many times while
the noise of the data, baseline, and integration bounds are varied based on a
noise model and user supplied probability distributions.

To view the documentation, click the badge below:

[![Documentation, latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://nluetts.github.io/NoisySignalIntegration.jl/dev/)

## Installation

The package is not yet registered in Julia's general package registry, you have to install it directly from this Github repository.

To install it for your project, enter the package mode in the Julia REPL (press `]`) and type:

```
add https://github.com/nluetts/NoisySignalIntegration.jl
```

While still in package mode, you can type

```
test NoisySignalIntegration
```

to run the package's unit tests.

## Getting Started

Check out the [documentation](https://nluetts.github.io/NoisySignalIntegration.jl/dev/) to learn how to use the package.

If you don't have a local Julia installation, you can test the package on [mybinder.org](https://mybinder.org).
Click the badge below and open one of the example notebooks:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nluetts/NSI-Binder/HEAD?filepath=example-1.ipynb&urlpath=lab)

## Contributing and Support

If you have problems with the package, submit an [issue](https://github.com/nluetts/NoisySignalIntegration.jl/issues).
Feel free to fork the project and open a pull request if you would like to contribute.