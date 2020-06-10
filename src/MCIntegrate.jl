module MCIntegrate

using Distributions: ContinuousUnivariateDistribution
using Distributions: MvNormal, Beta, LocationScale
using LinearAlgebra: eachcol
using LsqFit: curve_fit
using Plots
using Polynomials: fit
using Printf: @printf, @sprintf
using Random
using StatsBase: autocov

include("spectrum.jl")
include("noise.jl")
include("bounds.jl")
include("integration.jl")
include("mc.jl")

export Spectrum, crop, plot
export Noise, GaussianNoiseModel, MvGaussianNoiseModel, get_cov, estimate_autocov, fit_noise, plot_autocov, sample
export scale_shift_beta, LeftRightBound, WidthBound
export mc_integrate

# dev export
#export allapproxequal, detrend

end # module