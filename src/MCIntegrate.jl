module MCIntegrate

using Distributions: MvNormal, Beta
using LinearAlgebra: eachcol
using LsqFit: curve_fit
using Plots
using Polynomials: fit
using Printf: @printf, @sprintf
using Random: seed!
using StatsBase: autocov

include("spectrum.jl")
#include("analyze.jl")
include("noise.jl")
#include("integration.jl")
export Noise, GaussianNoiseModel, MvGaussianNoiseModel, Spectrum, crop, fit_noise, get_cov, plot, plot_autocov, sample, estimate_autocov

# dev export
export allapproxequal, detrend

end # module