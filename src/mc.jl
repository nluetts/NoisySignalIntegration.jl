function mc_integrate(
    s::Spectrum,
    nm::AbstractNoiseModel,
    bs::Vector{T}
    ;
    N::Int64=100_000
) where {T<:AbstractUncertainBound}
    integral_samples = Array{Float64}(undef, N, length(bs))

    noise_samples = sample(nm, length(s), N)
    bound_samples = [typeof(b) <: WidthBound ? sample(b, s + noise_samples[:, i], N) : sample(b, N) for (i, b) in enumerate(bs)]

    for i in 1:N
        spec = s + noise_samples[:, i]
        for j in 1:length(bs)
            left, right = bound_samples[j][i]
            integral_samples[i, j] = trapz(s.x, spec.y, left, right)
        end
    end

    return integral_samples
end # mc_integrate
