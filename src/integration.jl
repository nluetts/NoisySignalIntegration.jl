"""
    mc_integrate(uc::UncertainCurve{T, N}, bnds::Vector{UncertainBound{T, M}}; intfun=trapz)
    mc_integrate(uc::UncertainCurve{T, N}, bnds::UncertainBound{T, M}; intfun=trapz)

Integrate uncertain curve using uncertain bound(s).

Applies Monte-Carlo integration algorithm to samples stored in `uc` and returns one or an array
of uncertain areas of type `Particles{T, N}`.

# Keyword arguments

`intfun`:
The core integration function that is used to numerically integrate each draw.
Defaults to `NoisySignalIntegration.trapz`.
The function that is used to substitute [`trapz`](@ref) must share its call signature.

`subtract_baseline`:
If true, for each draw a local baseline defined by the integration window start and end point will be subtracted.
"""
function mc_integrate(uc::UncertainCurve{T, N}, bnds::Vector{UncertainBound{T, M}}; intfun=trapz, subtract_baseline=true) where {T, M, N}

    M != N && error("Samples sizes incompatible")
    
    areas = Array{T}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Integrating draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            l, r = get_draw(i, b)
            areas[i, j] = intfun(uc.x, cᵢ.y, l, r; subtract_baseline=subtract_baseline)
        end
    end
    return [Particles(areas[:,i]) for i in 1:size(areas)[2]]
end

function mc_integrate(uc::S, bnd::T; intfun=trapz, subtract_baseline=true) where {S <: UncertainCurve, T <: UncertainBound}
    return mc_integrate(uc, [bnd]; intfun=intfun, subtract_baseline=subtract_baseline)[1]
end