using MCIntegrate
using Random: seed!
using StatsBase: std
using Test

const mci = MCIntegrate

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


@testset "sample()" begin
    # sample GaussianNoiseModel
    let
        seed!(1)
        nm = GaussianNoiseModel(3.0)
        spls = nm |> sample
        @test size(spls) == (100,)
        @test size(sample(nm, 101, 2)) == (101, 2)
        @test std(sample(nm, 1_000_000)) ≈ 3.0 atol=1e-2
    end
    # sample MvGaussianNoiseModel
    let
        seed!(2)
        nm = MvGaussianNoiseModel(1.0, 3.0, 1.5)
        spls = nm |> sample
        @test size(spls) == (100,)
        @test size(sample(nm, 101, 2)) == (101, 2)
        # This step tests that the fit_noise function
        # actually reproduces the input parameters of
        # the noise generator model:
        noise_sample = sample(nm, 6000)[:,1]
        noise_param = noise_sample |> Noise |> fit_noise
        # The longer the noise sample is, the close the fit
        # results come to the input parameters (α = 3, λ = 1.5)
        # but increasing the noise sample length beyond 10k
        # points makes the test take _really_ long.
        # With 6k points (and the current seed value of 2),
        # the fitted parameters should reproduce the input
        # parameters with 5% accuray:
        @test noise_param.α ≈ 3.0 atol=5e-2
        @test noise_param.λ ≈ 1.5 atol=5e-2
    end
end