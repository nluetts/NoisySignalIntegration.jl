module MCIntegrate

using Distributions: ContinuousUnivariateDistribution
using Distributions: MvNormal, Beta, LocationScale, mean
using LinearAlgebra: eachcol
using LsqFit: curve_fit
using Plots
using Polynomials: fit
using Printf: @printf, @sprintf
using Random
using RecipesBase: @recipe
using StatsBase: autocov

include("common.jl")
include("types.jl")
include("stats.jl")
include("utils.jl")
include("plotting.jl")
include("mc.jl")

export
    Curve,
    Noise,
    GaussianNoiseModel,
    MvGaussianNoiseModel,
    LeftRightBound,
    WidthBound,
    scale_shift_beta,
    get_cov,
    estimate_autocov,
    fit_noise,
    plot_autocov,
    sample,
    sample!,
    crop,
    clone,
    mc_integrate

end # module