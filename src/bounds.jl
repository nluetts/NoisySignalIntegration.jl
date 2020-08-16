# Defines the UncertainBound type to handle uncertain untegration start and end points

struct UncertainBound{T, N}
    left::Particles{T, N}
    right::Particles{T, N}
end

# Constructors

"""
Create a left/right bound
"""
function UncertainBound(left::S, right::T, N::Int=10_000) where {S <: ContinuousUnivariateDistribution, T <: ContinuousUnivariateDistribution}
    left  = Particles(N, left)
    right = Particles(N, right)
    return UncertainBound(left, right)
end

"""
Create several, correlated width bounds
"""
function UncertainBound(
    pos::Vector{T},
    width::ContinuousUnivariateDistribution,
    uc::UncertainCurve{T, N}
) where {T, N}

    M = length(pos)
    left = Array{T}(undef, M, N)
    right = Array{T}(undef, M, N)
    width = Particles(N, width)

    for i ∈ 1:M
        pᵢ = pos[i]
        for j ∈ 1:N
            cⱼ = get_draw(j, uc) # sample j from curve
            wⱼ = get_draw(j, width) # sample j from width
            left[i, j], right[i, j] = left_right_from_peak(uc.x, cⱼ.y, pᵢ, wⱼ)
        end
    end
    
    return [UncertainBound(Particles(left[i, :]), Particles(right[i, :])) for i in 1:M]
end

"""
Create single width bound
"""
function UncertainBound(
    pos::T,
    width::ContinuousUnivariateDistribution,
    uc::UncertainCurve{T, N}
) where {T, N}
    bnd = UncertainBound([pos], width, uc)
    return bnd[1]
end


get_draw(n, bnd::UncertainBound) = [get_draw(n, bnd.left), get_draw(n, bnd.right)]


Statistics.mean(bnd::UncertainBound) = [Statistics.mean(bnd.left), Statistics.mean(bnd.right)]

"""
    scale_shift_beta(α, β, a, b)

Create a scaled and shifted `Beta(α, β)` distribution.
Samples fall in the interval [`a`, `b`].
"""
function scale_shift_beta(α, β, a, b)
    return LocationScale(a, b - a, Beta(α, β))
end
