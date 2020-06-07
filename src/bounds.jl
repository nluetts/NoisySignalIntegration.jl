
#######################################################
### Scaled and shifted distributions ##################
#######################################################

"""
    ScaledShiftedCUD

Scaled and shifted ContinuousUnivariateDistribution.

Note: Rescales and shifts random variable X to (X - offset) * scale + location
"""
struct ScaledShiftedCUD
    distribution::ContinuousUnivariateDistribution
    location::Float64
    scale::Float64
    offset::Float64
end

function sample(rng::AbstractRNG, d::ScaledShiftedCUD, dims::Integer)
    samples = rand(rng, d.distribution, dims)
    return (samples .- d.offset) .* d.scale .+ d.location
end

function sample(d::ScaledShiftedCUD, dims::Integer)
    return sample(Random.GLOBAL_RNG, d, dims)
end

#######################################################
### bounds datatype and sampling ######################
#######################################################

struct UncertainBounds{T<:ScaledShiftedCUD}
    left::T
    right::T
end

struct UncertainWidth{T<:ScaledShiftedCUD}
    width::T
end

function sample(b::UncertainBounds, samples::Integer)
    left = sample(b.left, samples)
    right = sample(b.right, samples)
    return hcat(left, right)
end
