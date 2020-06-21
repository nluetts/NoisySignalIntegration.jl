using Debugger

function mc_integrate(
    с::Curve,
    nm::AbstractNoiseModel,
    bs::Vector{T}
    ;
    N::Int64=100_000
) where {T<:AbstractUncertainBound}
    integral_samples = Array{Float64}(undef, N, length(bs))

    noise_samples = sample(nm, length(с), N)
    bound_samples = [typeof(b) <: WidthBound ? sample(b, с + noise_samples[:, i], N) : sample(b, N) for (i, b) in enumerate(bs)]

    for i in 1:N
        spec = с + noise_samples[:, i]
        for j in 1:length(bs)
            println(size(bound_samples[j]))
            left, right = bound_samples[j][i]
            integral_samples[i, j] = trapz(с.x, spec.y, left, right)
        end
    end

    return integral_samples
end # mc_integrate
