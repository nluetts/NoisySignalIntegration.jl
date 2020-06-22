function mc_integrate(
    crv::Curve,
    nm::AbstractNoiseModel,
    bs::Vector{T}
    ;
    N::Int64=100_000
) where {T<:AbstractUncertainBound}   
    # allocations
    m  = length(crv)
    nb = length(bs)
    integral_samples = Array{Float64}(undef, N, length(bs))
    noise_sample = Array{Float64}(undef, m)
    bound_sample = Array{Float64}(undef, 1, 2)

    for i in 1:N
        sample!(noise_sample, nm)
        spec = crv + noise_sample
        for (j, b) in enumerate(bs)
            typeof(b) == WidthBound ? sample!(bound_sample, b, c) : sample!(bound_sample, b)
            integral_samples[i, j] = trapz(crv.x, spec.y, bound_sample[1], bound_sample[2])
        end
    end

    return integral_samples
end # mc_integrate
