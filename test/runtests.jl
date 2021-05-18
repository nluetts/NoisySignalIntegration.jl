using Distributions
using NoisySignalIntegration
using MonteCarloMeasurements
using Random: seed!, rand
using Statistics
using StatsBase: mean, percentile, std
using Test

const NP = MonteCarloMeasurements.DEFAULT_NUM_PARTICLES
const nsi = NoisySignalIntegration

if false
    include("tmp.jl")
else
    include("test_common.jl")
    include("test_curves.jl")
    include("test_bounds.jl")
    include("test_noise.jl")
    include("test_integration.jl")
    include("test_regression.jl")
end