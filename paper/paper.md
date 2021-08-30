---
title: 'NoisySignalIntegration.jl: A Julia package for uncertainty evaluation of numeric integrals'
tags:
  - Julia
  - chemistry
  - physics
  - measurement uncertainty
  - noisy data
  - numeric integration
authors:
  - name: Nils O.&nbsp;B. Lüttschwager
    orcid: 0000-0001-8459-1714
    affiliation: 1
affiliations:
 - name: Georg-August-University Göttingen, Institute of Physical Chemistry, Tammannstraße 6, DE-37077 Göttingen
   index: 1
date: 8 July 2021
bibliography: paper.bib
---

# Summary

The evaluation of peak or band areas is a recurring task in scientific data
evaluation. For example, in molecular spectroscopy, absorption line or band
areas are often used to determine substance abundance.
NoisySignalIntegration.jl provides functionality to evaluate such signal
areas and associated uncertainties using a Monte-Carlo approach. Uncertainties
may include contributions from (potentially correlated) Gaussian noise,
baseline subtraction, and uncertainty in placing integration bounds. Uncertain
integration bounds can be defined in several ways to constrain the integration
based on the physical system under investigation (asymmetric signals, symmetric
signals, signals with identical width). The package thus offers a more
objective uncertainty evaluation than a statement based on experience or
laborious manual analysis [@goebench].

NoisySignalIntegration.jl includes a [detailed
documentation](https://nluetts.github.io/NoisySignalIntegration.jl/stable/) that
covers the typical workflow with several examples. The API uses custom
datatypes and convenience functions to aid the data
analysis and permits flexible customizations: Any probability distribution from
Distributions.jl [@distributions1; @distributions2] is a valid input to express uncertainty
in integration bounds, thus allowing to adapt the uncertainty analysis as
needed to ones state of knowledge. The core integration function can be swapped
if the included trapezoidal integration is deemed unsatisfactory in terms of
accuracy. The package uses MonteCarloMeasurements.jl [@mcm] to express
uncertain numbers which enables immediate uncertainty propagation.

# Statement of need

Several open source options for uncertainty propagation are available, e.g.
the Python packages uncertainties [@uncertainties] or MetroloPy [@metrolopy] or
Julia packages Measurements.jl [@measurements] and the aforementioned
MonteCarloMeasurements.jl [@mcm], but uncertainty evaluation when integrating
experimental x-y data is not fully addressed and requires significant effort
from the user. A straightforward solution to this problem is to fit a peak
function of appropriate shape and derive the uncertainty from the fit. However,
the data may not be described by a peak function and/or the noise may not be
normally distributed, preventing a simple least squares regression.
NoisySignalIntegration.jl was developed as a solution to this problem.  While
the package was developed specifically for the determination of band area
uncertainties in the context of molecular spectroscopy [@karir2019; @gawrilow; @zimmermann],
it is applicable in any research area where signals (peaks, lines, bands,
etc.)^[The name usually depends on the specific area and context.] in x-y data
need to be integrated and a thorough uncertainty analysis is desired.

# Acknowledgements

The author thanks Charlotte Zimmermann, Maxim Gawrilow, and Robert Medel for testing and
discussions during the development of NoisySignalIntegration.jl.

# References
