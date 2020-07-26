# MCIntegrate.jl

*MCIntegralte.jl -- A tool to determine uncertainty of numeric integrals of noisy data.*

`MCIntegrate` is a package that implements a method to determine the uncertainty of numeric integrals of noisy data on the basis of a Monte-Carlo process.
The package was originally intended to estimate uncertainties of band signals in FTIR
spectra, which is reflected in the examples given in the following usage guide.



```@docs
get_cov
crop
```