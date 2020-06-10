function mc_integrate(s::Spectrum, nm::NoiseModel, bs::Vector{AbstractUncertainBound}; N::Int64=100_000)

    integral_samples = Array{Float64}(undef, N, length(bs))

    for i in 1:N
        for (j, bound) in enumerate(bs)
            spec = s + sample(nm, length(s))
            left, right = if typeof(bound) == UncertainWidth
                    sample(bound, spec)
                else
                    sample(bound)
                end
            integral_samples[i, j] = trapz(s.x, spec.y, left, right)
        end
    end

    return integral_samples
end # mc_integrate
