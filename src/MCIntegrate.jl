module MCIntegrate

using Distributions:        ContinuousUnivariateDistribution
using Distributions:        MvNormal, Beta, LocationScale, mean
using LinearAlgebra:        eachcol
using LsqFit:               curve_fit
using MonteCarloMeasurements: Particles
using Plots
using Polynomials:          fit
using Printf:               @printf, @sprintf
using Random
using RecipesBase:          @recipe
using StatsBase:            autocov

include("common.jl")
include("curves.jl")
include("bounds.jl")
include("noise.jl")
include("plotting.jl")
include("integration.jl")

export
    Curve,
    GaussianNoiseModel,
    MvGaussianNoiseModel,
    NoiseSample,
    UncertainCurve,
    UncertainBound,
    add_noise,
    correlated_noise,
    crop,
    fit_noise,
    get_cov,
    mc_integrate,
    plot,
    plot_autocov,
    scale_shift_beta,
    uncorrelated_noise

end # module