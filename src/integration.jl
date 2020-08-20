"""
    mc_integrate(uc::UncertainCurve{T, N}, bnds::Vector{UncertainBound{T, M}}; intfun=trapz)
    mc_integrate(uc::UncertainCurve{T, N}, bnds::UncertainBound{T, M}; intfun=trapz)

Integrate uncertain curve using uncertain bound(s).

Applies Monte-Carlo integration algorithm to samples stored in `uc` and returns one or an array
of uncertain areas of type `Particles{T, N}`.
"""

function mc_integrate(uc::UncertainCurve{T, N}, bnds::Vector{UncertainBound{T, M}}; intfun=trapz) where {T, M, N}

    M != N && error("Samples sizes incompatible")
    
    areas = Array{T}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Integrating draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            l, r = get_draw(i, b)
            areas[i, j] = intfun(uc.x, cᵢ.y, l, r)
        end
    end
    return [Particles(areas[:,i]) for i in 1:size(areas)[2]]
end

mc_integrate(uc::S, bnd::T; intfun=trapz) where {S <: UncertainCurve, T <: UncertainBound} = mc_integrate(uc, [bnd]; intfun=intfun)[1]