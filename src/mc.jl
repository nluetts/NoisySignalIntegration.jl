function mc_integrate(us::UncertainCurve{T, N}, bnds::Vector{UncertainBound{T, M}}) where {T, M, N}
    M != N && error("Samples sizes incompatible")
    
    areas = Array{Float64}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 100 == 0 && print("Integrating draw $i/$N \r")
        y = get_draw(i, us)
        for (j, b) in enumerate(bnds)
            l, r = get_draw(i, b)
            areas[i, j] = MCIntegrate.trapz(us.x, y, l, r)
        end
    end
    return [Particles(areas[:,i]) for i in 1:size(areas)[2]]
end