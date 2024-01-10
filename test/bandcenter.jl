module BandCenter

using Distributions
using Random: seed!
using MonteCarloMeasurements
using NoisySignalIntegration
using NoisySignalIntegration: lininterp, _local_baseline, _endpoint_to_endpoint_baseline, get_draw
using Plots
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

function baseline_from_points(x1, x2, y1, y2)
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    @assert b ≈ y2 - m * x2 rtol=1e-10
    m, b
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
    m, b = baseline_from_points(left, right, yleft, yright)
    band_center(curve, left, right, (m, b))
end


function mc_bandcenter(uc::UncertainCurve{T,N}, bnds::Vector{UncertainBound{T,M}}; subtract_baseline=false, local_baseline=false) where {T,M,N}

    M != N && error("Samples sizes incompatible")
    subtract_baseline && @warn("subtract_baseline keyword argument is deprecated, use local_baseline instead.")
    (subtract_baseline && local_baseline) && error("local_baseline and subtract_baseline cannot both be true.") |> throw

    centers = Array{T}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Processing draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            xₗ, xᵣ = get_draw(i, b)
            if local_baseline
                baseline = baseline_from_points(_local_baseline(cᵢ.x, cᵢ.y, xₗ, xᵣ, b)...)
                centers[i, j] = band_center(cᵢ, xₗ, xᵣ, baseline)
            elseif subtract_baseline
                centers[i, j] = band_center(cᵢ, xₗ, xᵣ, true)
            else
                centers[i, j] = band_center(cᵢ, xₗ, xᵣ)

            end
        end
    end
    return [Particles(centers[:, i]) for i in 1:size(centers)[2]]
end

function mc_bandcenter(uc::S, bnd::T; subtract_baseline=false, local_baseline=false) where {S<:UncertainCurve,T<:UncertainBound}
    return mc_bandcenter(uc, [bnd]; subtract_baseline=subtract_baseline, local_baseline=local_baseline)[1]
end

function main()

    spectrum = NoisySignalIntegration.testdata_1()
    slice_bands = crop(spectrum,  5.0,  40.0)
    slice_noise = crop(spectrum, 40.0, 100.0)

    plot(slice_bands; label="bands")
    plot!(slice_noise; label="noise")

    noise = NoiseSample(slice_noise, 3)
    nm = fit_noise(noise)
    uncertain_spectrum = add_noise(slice_bands, nm)
    position = [15.0, 30.0]
     # widths will fall in the range 2 to 3, with a maximum at 2.5
    width_distribution = scale_shift_beta(2, 2, 3, 4)
    # define a "width bound"
    bds = UncertainBound(position, width_distribution, uncertain_spectrum)
    ctrs = mc_bandcenter(uncertain_spectrum, bds, local_baseline=true)
    mc_bandcenter(uncertain_spectrum, bds, subtract_baseline=true)
    mc_bandcenter(uncertain_spectrum, bds)
    uncertain_spectrum, bds, ctrs
end

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
