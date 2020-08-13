using Test
using MCIntegrate
using Random: seed!

const mci = MCIntegrate

@testset "detrend()" begin
    for n ∈ 2:4
        @test begin
            x = collect(Float64, 1:20)
            y = x.^n
            all(yᵢ + 1.0 ≈ 1.0 for yᵢ ∈ mci.detrend(x, y, n)) == true
        end
    end
    @test begin
        x = collect(Float64, 1:20)
        y = @. 0.1 + 2x - 0.01x^2
        all(yᵢ + 1.0 ≈ 1.0 for yᵢ ∈ mci.detrend(x, y, 2)) == true
    end
end


@testset "gauss_kernel()" begin
    @test mci.gauss_kernel(0, [2, 1]) == 4.0
    @test mci.gauss_kernel(2.0, [2, 2]) == 4.0 * exp(-0.5)
    @test mci.gauss_kernel(4.0, [2, 1]) == 4.0 * exp(-8)
    @test mci.gauss_kernel(4.0, [0, 1]) == 0.0
end


@testset "get_cov()" begin
    let f = 1.0000001
        @test_throws ArgumentError mci.get_cov(1.0, -1, 2.0, 1.0)
        @test mci.get_cov(1.0, 3, 2.0, 1.0) == 4 .* [f         exp(-0.5) exp(-2)  ;
                                                     exp(-0.5) f         exp(-0.5);
                                                     exp(-2)   exp(-0.5) f         ]
    end
end


@testset "correlated_noise()" begin
    seed!(42)
    nm = MvGaussianNoiseModel(1.0, 1.0, 5.0)
    cnoise = correlated_noise(nm, 1000, 3)
    @test length(cnoise) == 1000
    @test length(cnoise[1].particles) == 3
    cnoise_sample = NoiseSample(mean(cnoise))
    fit_noise(cnoise_sample) |> println
end

@testset "NoiseSample" begin

    # test basic operation

    noise = NoiseSample([1.0, 2.0], [4.0, 4.0])
    @test length(noise) == 2
    @test noise[end][1] == 2.0
    @test noise[end][2] + 1 ≈ 1.0
    @test [noise; noise] == NoiseSample([1.0, 2.0, 1.0, 2.0], [4.0, 4.0, 4.0, 4.0])

    # x and y need to have same length
    @test_throws ArgumentError NoiseSample([1.0], [4.0, 4.0])
    @test_throws ArgumentError NoiseSample([1.0, 2.0], [4.0])
end;