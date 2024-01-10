module BandCenter

using NoisySignalIntegration
using NoisySignalIntegration: lininterp
using Test


function test_data_1()
    # gaussian curve generated in LibreOffice
    X = [
        0	4.53999297624849E-05;
        1	0.000303539138078867;
        2	0.00166155727317393;
        3	0.00744658307092434;
        4	0.0273237224472926;
        5	0.0820849986238988;
        6	0.201896517994655;
        7	0.406569659740599;
        8	0.670320046035639;
        9	0.90483741803596;
        10	1;
        11	0.90483741803596;
        12	0.670320046035639;
        13	0.406569659740599;
        14	0.201896517994655;
        15	0.0820849986238988;
        16	0.0273237224472926;
        17	0.00744658307092434;
        18	0.00166155727317393;
        19	0.000303539138078867;
        20	4.53999297624849E-05;
    ]
    Curve(X[:, 1], X[:, 2])
end

function band_center(
    curve::Curve{T},
    left::T,
    right::T,
    baseline::Union{Nothing, Tuple{T, T}}=nothing
) where {T<:AbstractFloat}
    sum_intensity = zero(T)
    sum_weighted_position = zero(T)
    m, b = isnothing(baseline) ? (zero(T), zero(T)) : baseline
    for (xi, yi) in (curve, left, right)
        sum_intensity += yi - m * xi - b
        sum_weighted_position += xi * (yi - m * xi - b)
    end
    sum_weighted_position / sum_intensity
end

# this will behave like `local_baseline`, i.e. spanning a linear
# baseline from the start to the endpoint
function band_center(
    curve::Curve{T},
    left::T,
    right::T,
    baseline::Bool
) where {T<:AbstractFloat}
    yleft = lininterp(left, curve)
    yright = lininterp(right, curve)
    m = (yright - yleft) / (right - left)
    b = yleft - m * left
    @assert b == yright - m * right
    band_center(curve, left, right, (m, b))
end


# function mc_bandcenter(uc::UncertainCurve{T,N}, bnds::Vector{UncertainBound{T,M}}; intfun=trapz, subtract_baseline=false, local_baseline=false) where {T,M,N}

#     M != N && error("Samples sizes incompatible")
#     subtract_baseline && @warn("subtract_baseline keyword argument is deprecated, use local_baseline instead.")
#     (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

#     areas = Array{T}(undef, N, length(bnds))
#     for i ∈ 1:N
#         i % 1000 == 0 && print("Integrating draw $i/$N \r")
#         cᵢ = get_draw(i, uc)
#         for (j, b) in enumerate(bnds)
#             xₗ, xᵣ = get_draw(i, b)
#             xs, ys = uc.x, cᵢ.y
#             areas[i, j] = intfun(xs, ys, xₗ, xᵣ)
#             if local_baseline
#                 areas[i, j] -= singletrapz(_local_baseline(xs, ys, xₗ, xᵣ, b)...)
#             end
#             if subtract_baseline
#                 areas[i, j] -= singletrapz(_endpoint_to_endpoint_baseline(xs, ys, xₗ, xᵣ)...)
#             end
#         end
#     end
#     return [Particles(areas[:, i]) for i in 1:size(areas)[2]]
# end

# function mc_bandcenter(uc::S, bnd::T; intfun=trapz, subtract_baseline=false, local_baseline=false) where {S<:UncertainCurve,T<:UncertainBound}
#     return mc_bandcenter(uc, [bnd]; intfun=intfun, subtract_baseline=subtract_baseline, local_baseline=local_baseline)[1]
# end

function run_tests()
    crv = test_data_1()
    @testset "gaussian peak" begin
        # without baseline
        # @test band_center(crv, 0.0, 20.0) == 10.0
        # @test band_center(crv, 6.5, 13.5) == 10.0
        # @test band_center(crv, 6.5, 15.5) ≈ 10.0809916179 atol=1e-10
        # @test band_center(crv, 6.3, 11.85) ≈ 9.540840963609 atol=1e-10
        # with baseline subtraction
        # @test band_center(crv, 6.3, 11.85, (0.025, -0.43)) ≈ 9.3215429658485 atol=1e-10
        @test band_center(crv, 6.3, 11.85, true) ≈ 9.3215429658485 atol=1e-10
    end
end

end
