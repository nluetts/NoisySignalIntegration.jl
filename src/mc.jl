function mc_integrate(
    s::Spectrum,
    nm::AbstractNoiseModel,
    bs::Vector{T}
    ;
    N::Int64=100_000
) where {T<:AbstractUncertainBound}
    integral_samples = Array{Float64}(undef, N, length(bs))

    spectrum_samples = s.y .+ sample(nm, length(s), N)

    for i in 1:N
        #println(i)
        spec = spectrum_samples[:, i]
        #println(spec)
        for (j, bound) in enumerate(bs)
            left, right = if typeof(bound) <: WidthBound
                    sample(bound, Spectrum(s.x, spec))
#                    println(spl)
#                    spl
                else
#                    println(typeof(bound))
                    sample(bound)
#                    println(spl)
#                    spl
                end
            integral_samples[i, j] = trapz(s.x, spec, left, right)
        end
    end

    return integral_samples
end # mc_integrate
