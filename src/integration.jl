function _local_baseline(x, y, l, r, b::UncertainBound)
    ŷₗ = mean(lininterp(q, x, y) for q in b._left_quantiles)
    ŷᵣ = mean(lininterp(q, x, y) for q in b._right_quantiles)
    x̂ₗ = b._left_quantiles[IDX_BOUND_MEDIAN]
    x̂ᵣ = b._right_quantiles[IDX_BOUND_MEDIAN]
    yl = lininterp(l, x̂ₗ, x̂ᵣ, ŷₗ, ŷᵣ)
    yr = lininterp(r, x̂ₗ, x̂ᵣ, ŷₗ, ŷᵣ)
    return l, r, yl, yr
end


_endpoint_to_endpoint_baseline(x, y, l, r) = (l, r, lininterp(l, x, y), lininterp(r, x, y))


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
function mc_integrate(uc::UncertainCurve{T, N}, bnds::Vector{UncertainBound{T, M}}; intfun=trapz, subtract_baseline=false, local_baseline=false) where {T, M, N}

    M != N && error("Samples sizes incompatible")
    println(local_baseline, subtract_baseline)
    subtract_baseline && @warn("subtract_baseline keyword argument is deprecated, use local_baseline instead.")
    (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

    areas = Array{T}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Integrating draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            l, r = get_draw(i, b)
            x, y = uc.x, cᵢ.y
            areas[i, j] = intfun(x, y, l, r)
            if local_baseline
                areas[i, j] -= singletrapz(_local_baseline(x, y, l, r, b)...)
            end
            if subtract_baseline
                areas[i, j] -= singletrapz(_endpoint_to_endpoint_baseline(x, y, l, r)...)
            end
        end
    end
    return [Particles(areas[:,i]) for i in 1:size(areas)[2]]
end

function mc_integrate(uc::S, bnd::T; intfun=trapz, subtract_baseline=false, local_baseline=false) where {S <: UncertainCurve, T <: UncertainBound}
    return mc_integrate(uc, [bnd]; intfun=intfun, subtract_baseline=subtract_baseline, local_baseline=local_baseline)[1]
end