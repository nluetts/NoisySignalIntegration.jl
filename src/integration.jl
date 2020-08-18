function mc_integrate(
    us::UncertainCurve{T, N},
    bnds::Vector{UncertainBound{T, M}}
    ;
    intfun=trapz
) where {T, M, N}

    M != N && error("Samples sizes incompatible")
    
    areas = Array{T}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Integrating draw $i/$N \r")
        cᵢ = get_draw(i, us)
        for (j, b) in enumerate(bnds)
            l, r = get_draw(i, b)
            areas[i, j] = intfun(us.x, cᵢ.y, l, r)
        end
    end
    return [Particles(areas[:,i]) for i in 1:size(areas)[2]]
end