module IterAllocTest
using NoisySignalIntegration
using NoisySignalIntegration: get_draw, _local_baseline, _endpoint_to_endpoint_baseline, singletrapz
using MonteCarloMeasurements
using Random: seed!
using Test

function integrate_draw(crv::Curve{T}, left::T, right::T) where {T}
    area = zero(T)
    for ((xi, yi), (xj, yj)) in zip((crv, left, right), Iterators.drop((crv, left, right), 1))
        area += singletrapz(xi, xj, yi, yj)
    end
    area
end

function mc_integrate_2(uc::UncertainCurve{T,N}, bnds::Vector{UncertainBound{T,M}}; intfun=trapz, subtract_baseline=false, local_baseline=false) where {T,M,N}

    M != N && error("Samples sizes incompatible")
    subtract_baseline && @warn("subtract_baseline keyword argument is deprecated, use local_baseline instead.")
    (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

    areas = Array{T}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Integrating draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            xₗ, xᵣ = get_draw(i, b)
            areas[i, j] = integrate_draw(cᵢ, xₗ, xᵣ)
            if local_baseline
                areas[i, j] -= singletrapz(_local_baseline(cᵢ.x, cᵢ.y, xₗ, xᵣ, b)...)
            end
            if subtract_baseline
                areas[i, j] -= singletrapz(_endpoint_to_endpoint_baseline(cᵢ.x, cᵢ.y, xₗ, xᵣ)...)
            end
        end
    end
    return [Particles(areas[:, i]) for i in 1:size(areas)[2]]
end

function mc_integrate_2(uc::S, bnd::T; intfun=trapz, subtract_baseline=false, local_baseline=false) where {S<:UncertainCurve,T<:UncertainBound}
    return mc_integrate_2(uc, [bnd]; intfun=intfun, subtract_baseline=subtract_baseline, local_baseline=local_baseline)[1]
end

function get_test_spectrum(seed)
    seed!(seed)
    x = collect(0:0.1:200)
    y = @. exp(-(x - 15)^2) + 2 * exp(-(x - 30)^2)
    @. y += 1.0 + x * 1.5e-2 - (x - 50)^3 * 3e-7 # add polynomial baseline
    n = length(x)
    δx = x[2] - x[1]
    return Curve(x, y .+ (get_cov(δx, n, 0.1, 0.5) |> MvNormal |> rand))
end

function main()
    c = crop(get_test_spectrum(1), 10, 40)
    uc = add_noise(c, MvGaussianNoiseModel(0.1, 0.1, 0.5))
    ub = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 3.0, 4.0), uc)
    @time area = mc_integrate_2(uc, ub, local_baseline=true)
    @time area = mc_integrate(uc, ub, local_baseline=true)
end


end
