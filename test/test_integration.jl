function get_test_spectrum(seed)
    seed!(seed)
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
        area = mc_integrate(uc, ub, subtract_baseline=true)
        @test mean(area) ≈ 1.808961690 atol = 1e-8
        @test std(area) ≈ 0.2488841198 atol = 1e-8
    end
    @testset "integration of two peaks" begin
        c = crop(get_test_spectrum(1), 10, 40)
        uc = add_noise(c, MvGaussianNoiseModel(0.1, 0.1, 0.5))
        ubs = UncertainBound([15., 30.], scale_shift_beta(2.0, 2.0, 3.0, 4.0), uc)
        area1, area2 = mc_integrate(uc, ubs, subtract_baseline=true)
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

@testset "integrate using Simpsons rule from external package" begin

    using NumericalIntegration: integrate, SimpsonEven
    
    function simps(x, y, a, b; subtract_baseline=true)
        idx = [(xᵢ >= a) & (xᵢ <= b) for xᵢ in x]
        return integrate(x[idx], y[idx], SimpsonEven())
    end
    c = crop(get_test_spectrum(1), 10, 40)
    uc = add_noise(c, MvGaussianNoiseModel(0.1, 0.1, 0.5))
    ub = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 3.0, 4.0), uc)
    area_s = mc_integrate(uc, ub; intfun=simps)
    @test mean(area_s) ≈ 5.6868358 atol = 1e-7
    @test std(area_s) ≈ 0.32897085 atol = 1e-7
end

@testset "error if local_baseline and subtract_baseline true" begin
    c = NoisySignalIntegration.testdata_1()
    uc = add_noise(c, MvGaussianNoiseModel(0.1, 0.1, 0.5))
    ub = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 3.0, 4.0), uc)
    @test_throws ErrorException area = mc_integrate(uc, ub; local_baseline=true, subtract_baseline=true)
end


@testset "test _local_baseline" begin
    xs = 1.0:10.0 |> collect
    b = UncertainBound(Uniform(1.5, 2.5), Uniform(4.5, 6.5))
    @testset "test _local_baseline with linear dataset" begin
        ys = xs
        xₗ, xᵣ, yₗ, yᵣ = nsi._local_baseline(xs, ys, 1.75, 5.25, b)
        yyₗ = nsi.lininterp(1.75, xs, ys)
        yyᵣ = nsi.lininterp(5.25, xs, ys)
        @test yyₗ == yₗ
        @test yyᵣ == yᵣ
    end
    @testset "test _local_baseline with parabolic dataset" begin
        ys = xs.^2
        xₗ, xᵣ, yₗ, yᵣ = nsi._local_baseline(xs, ys, 2.25, 5.25, b)
        @test yₗ ≈ 6.1329 atol=1e-4  # results derived by doing computation by hand in Numpy
        @test yᵣ ≈ 28.8471 atol=1e-4
    end
    @testset "test _local_baseline with normal and uniform distribution on noise" begin
        noise = [
            1 1.51
            2 1.41
            3 1.2
            4 0.88
            5 0.5
            6 0.09
            7 -0.32
            8 -0.67
            9 -0.94
            10 -1.09
            11 -1.12
            12 -1.04
            13 -0.86
            14 -0.63
            15 -0.39
            16 -0.18
            17 -0.04
            18 0.02
            19 -0.01
            20 -0.12
            21 -0.25
            22 -0.36
            23 -0.43
            24 -0.43
            25 -0.4
            26 -0.36
            27 -0.35
            28 -0.39
            29 -0.48
            30 -0.6
            31 -0.73
            32 -0.8
            33 -0.8
            34 -0.71
            35 -0.53
            36 -0.3
            37 -0.06
            38 0.15
            39 0.27
            40 0.29
        ]
        bnd = UncertainBound(Normal(4.0, 2.0), Uniform(30.0, 35.0))
        xₗ, xᵣ, yₗ, yᵣ = nsi._local_baseline(noise[:,1], noise[:,2], 4.0, 32.5, bnd)
        @test yₗ ≈ 0.77916842211701 # manually derived in LibreCalc
        @test yᵣ ≈ -0.7210156
    end
end