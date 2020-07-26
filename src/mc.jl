WidthBoundUnion = Union{WidthBound, WidthBoundClone}

function _mc_integrate!(integral_samples, curve::Curve, noise_samples, bounds, N)
    for i in 1:N
        i % 100 == 0 && print("Integrating draw $i/$N \r")
        curve_ = curve + noise_samples[:, i]
        for (j, b) in enumerate(bounds)
            bound_sample = typeof(b) <: WidthBoundUnion ? sample(b, curve) : sample(b)
            integral_samples[i, j] = trapz(curve.x, curve_.y, bound_sample[1], bound_sample[2])
        end
    end
    println()
end


function mc_integrate(
    crv::Curve,
    bs::Vector{T},
    nm::AbstractNoiseModel
    ;
    N::Int64=100_000
) where {T<:AbstractUncertainBound}   

    integral_samples = Array{Float64}(undef, N, length(bs))

    # sample noise and bounds
    print("\nPreparing noise samples ... ")
    noise_samples = sample(nm, length(crv), N)
    print("done.\n")
    
    _mc_integrate!(integral_samples, crv, noise_samples, bs, N)

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
    
    _mc_integrate!(integral_samples, crv, noise_samples, bs, N)

    return integral_samples
end # mc_integrate
