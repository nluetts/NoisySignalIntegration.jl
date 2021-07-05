module NoisySignalIntegration

using Distributions:        ContinuousUnivariateDistribution
using Distributions:        MvNormal, Beta, LocationScale
using LsqFit:               curve_fit
using MonteCarloMeasurements
using Polynomials:          fit
using Printf:               @sprintf
using Random:               seed!
using RecipesBase
using Requires
using Statistics:           mean
using StatsBase:            autocov

import Statistics


include("common.jl")
include("curves.jl")
include("bounds.jl")
include("noise.jl")
include("integration.jl")
include("plotting.jl")
include("testdata.jl")

function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("animate.jl")
end

export
    Curve,
    GaussianNoiseModel,
    MvGaussianNoiseModel,
    NoiseSample,
    UncertainBound,
    UncertainCurve,
    add_noise,
    crop,
    stitch,
    fit_noise,
    get_cov,
    mc_integrate,
    mean,
    plotautocovfit,
    plot,
    scale_shift_beta,
    std,
    trapz,
    @samples

end # module