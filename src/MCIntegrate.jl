module MCIntegrate

using Distributions:        ContinuousUnivariateDistribution
using Distributions:        MvNormal, Beta, LocationScale, mean
using LinearAlgebra:        eachcol
using LsqFit:               curve_fit
using MonteCarloMeasurements
using Plots
using Polynomials:          fit
using Printf:               @printf, @sprintf
using Random:               seed!
using RecipesBase:          @recipe
using StatsBase:            autocov

import Statistics


include("common.jl")
include("curves.jl")
include("bounds.jl")
include("noise.jl")
include("plotting.jl")
include("integration.jl")
include("testdata.jl")

export
    Curve,
    GaussianNoiseModel,
    MvGaussianNoiseModel,
    NoiseSample,
    UncertainBound,
    UncertainCurve,
    add_noise,
    crop,
    fit_noise,
    get_cov,
    mc_integrate,
    mean,
    plot_autocov,
    plot,
    scale_shift_beta,
    std,
    trapz

end # module