#######################################################
### Scaled and shifted Beta ###########################
#######################################################

"""
    ScaledShiftedBeta

Scaled and shifted `Beta(α, β)` distribution.
Samples fall in the interval [`a`, `b`].
"""
ScaledShiftedBeta = LocationScale{Float64, Beta{Float64}}

"""
    scale_shift_beta(α::Float64, β::Float64, a::Float64, b::Float64)

Create a scaled and shifted `Beta(α, β)` distribution.
Samples fall in the interval [`a`, `b`].
"""
function scale_shift_beta(α::Float64, β::Float64, a::Float64, b::Float64)
    return LocationScale(a, b - a, Beta(α, β))
end

#######################################################
### bounds datatype and sampling ######################
#######################################################

abstract type AbstractUncertainBound end

struct LeftRightBound{S1<:ContinuousUnivariateDistribution,
                      S2<:ContinuousUnivariateDistribution} <: AbstractUncertainBound
    left::S1
    right::S2
end

struct WidthBound{T<:ContinuousUnivariateDistribution} <: AbstractUncertainBound
    loc::AbstractFloat
    width::T
end

function sample(b::LeftRightBound, samples::Integer)
    left = rand(b.left, samples)
    right = rand(b.right, samples)
    return hcat(left, right)
end


"""
    left_right_from_peak(x, y, p, w)

Find peak in interval `p ± w/2` and return symmetric bounds
left and right from peak with width `w`.
"""
function left_right_from_peak(x, y, p, w)
    # find indices of interval [left, right]
    left = p - w/2
    right = p + w/2
    l = searchsortedlast(x, left)    # left index:  x[l] <= p - w
    r = searchsortedfirst(x, right)  # right index: x[r] >= p + w
    
    m = l + argmax(view(y, l:r)) - 1 # index of local maximum at x[m]
    return [x[m] - w/2, x[m] + w/2]
end

sample(w::WidthBound, samples::Integer=1) = rand(w.width, samples)
function sample(w::WidthBound, s::Spectrum, samples::Integer=1) 
    ws = rand(w.width, samples)
    return [left_right_from_peak(s.x, s.y, w.loc, w) for w in ws]
end
    