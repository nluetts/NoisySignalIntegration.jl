module MCIntegrate

using Distributions: MvNormal, Beta, ContinuousUnivariateDistribution
using LinearAlgebra: eachcol
using LsqFit: curve_fit
using Plots
using Polynomials: fit
using Printf: @printf, @sprintf
using Random
using StatsBase: autocov

#include("spectrum.jl")
#include("noise.jl")
include("bounds.jl")
#export Spectrum, crop, plot
#export Noise, GaussianNoiseModel, MvGaussianNoiseModel, fit_noise, get_cov, plot_autocov, sample, estimate_autocov
export ScaledShiftedCUD, UncertainBounds, sample

# dev export
#export allapproxequal, detrend

end # module