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
    GaussianNoiseModel,
    LeftRightBound,
    MvGaussianNoiseModel,
    Noise,
    WidthBound,
    clone,
    crop,
    get_cov,
    mc_integrate,
    plot_autocov,
    sample,
    sample!,
    scale_shift_beta

end # module