using Test
using Distributions
using MCIntegrate
using Random


function get_test_spectrum(seed)
    Random.seed!(seed)
    x = collect(0:0.1:200)
    y = @. exp(-(x-15)^2) + 2 * exp(-(x-30)^2)
    @. y += 1.0 + x*1.5e-2 - (x-50)^3*3e-7 # add polynomial baseline
    n = length(x)
    δx = x[2] - x[1]
    return Curve(x, y .+ (get_cov(δx, n, 0.1, 0.5) |> MvNormal |> rand))
end

@testset "regression tests" begin
    @testset "integration of single peak" begin
        c = crop(get_test_spectrum(1), 10, 40)
        uc = add_noise(c, MvGaussianNoiseModel(0.1, 0.1, 0.5))
        ub = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 3.0, 4.0), uc)
        area = mc_integrate(uc, ub)
        @test mean(area) ≈ 1.808961690 atol = 1e-8
        @test std(area) ≈ 0.2488841198 atol = 1e-8
    end
    @testset "integration of two peaks" begin
        c = crop(get_test_spectrum(1), 10, 40)
        uc = add_noise(c, MvGaussianNoiseModel(0.1, 0.1, 0.5))
        ubs = UncertainBound([15., 30.], scale_shift_beta(2.0, 2.0, 3.0, 4.0), uc)
        area1, area2 = mc_integrate(uc, ubs)
        @test mean(area1) ≈ 1.808961690 atol = 1e-8
        @test std(area1) ≈ 0.2488841198 atol = 1e-8
        @test mean(area2) ≈ 3.442573414 atol = 1e-8
        @test std(area2) ≈ 0.3122904192 atol = 1e-8
    end
end

@testset "pass custom integration function" begin
    f(args...; kw...) = 0.0 # this should make integration yield area = 0 ± 0
    c = crop(get_test_spectrum(1), 10, 40)
    uc = add_noise(c, MvGaussianNoiseModel(0.1, 0.1, 0.5))
    ub = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 3.0, 4.0), uc)
    area = mc_integrate(uc, ub; intfun=f)
    @test mean(area) == 0.0
    @test std(area) == 0.0
end