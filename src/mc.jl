function mc_integrate(
    crv::Curve,
    nm::AbstractNoiseModel,
    bs::Vector{T}
    ;
    N::Int64=100_000
) where {T<:AbstractUncertainBound}   

    integral_samples = Array{Float64}(undef, N, length(bs))
    bound_sample = Array{Float64}(undef, 1, 2)

    # sample noise and bounds
    print("\nPreparing noise samples ... ")
    noise_samples = sample(nm, length(crv), N)
    print("done.\n")
    
    for i in 1:N
        i % 100 == 0 && print("Integrating draw $i/$N \r")
        spec = crv + noise_samples[:, i]
        for (j, b) in enumerate(bs)
            typeof(b) <: WidthBound ? sample!(bound_sample, b, crv) : sample!(bound_sample, b)
            integral_samples[i, j] = trapz(crv.x, spec.y, bound_sample[1], bound_sample[2])
        end
    end
    println()

    return integral_samples
end # mc_integrate


function mc_integrate(
    crv::Curve,
    noise_samples::Array{Float64, 2},
    bs::Vector{T}
    ;
    N::Int64=100_000
) where {T<:AbstractUncertainBound}   

    integral_samples = Array{Float64}(undef, N, length(bs))
    bound_sample = Array{Float64}(undef, 1, 2)
    
    for i in 1:N
        i % 100 == 0 && print("Integrating draw $i/$N \r")
        spec = crv + noise_samples[:, i]
        for (j, b) in enumerate(bs)
            typeof(b) <: WidthBound ? sample!(bound_sample, b, crv) : sample!(bound_sample, b)
            integral_samples[i, j] = trapz(crv.x, spec.y, bound_sample[1], bound_sample[2])
        end
    end
    println()

    return integral_samples
end # mc_integrate
