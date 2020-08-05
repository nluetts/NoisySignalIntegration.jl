using Test
using MCIntegrate

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
        spls = GaussianNoiseModel(3.0) |> sample
        @test size(spls) == (100, 1)
    end
end