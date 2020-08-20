using Test
using MCIntegrate
using Random: seed!
using StatsBase: mean

const mci = MCIntegrate

@testset "detrend()" begin
    # use detrend on polynomials and see if results are close to 0
    # (add 1. and see if results ≈ 1. as a workaround)
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


@testset "NoiseSample" begin

    # test construction

    let
        ns1 = NoiseSample{Int}([1, 2, 3], [4, 5, 6])
        ns2 = NoiseSample([1, 2, 3], [4, 5, 6])
        ns3 = NoiseSample(Curve([1, 2, 3], [4, 5, 6]))
        ns4 = NoiseSample([4, 5, 6])
        for ns in (ns1, ns2, ns3, ns4)
            @test length(ns) == 3
            @test eltype(ns) == Int
            @test ns.x == [1, 2, 3]
            @test ns.y == [-1, 0, 1]
        end
    end

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


@testset "GaussianNoiseModel" begin
    @test eltype(GaussianNoiseModel(1)) == Int
    @test eltype(GaussianNoiseModel(1.)) == Float64
    @test_throws MethodError GaussianNoiseModel(1 + 1im)
end


@testset "MvGaussianNoiseModel" begin
    @test eltype(MvGaussianNoiseModel(1, 1, 1)) == Int
    @test eltype(MvGaussianNoiseModel(1., 1., 1.)) == Float64
    @test_throws MethodError MvGaussianNoiseModel(1 + 1im, 1, 1)
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


@testset "add_noise()" begin
    let # use MvGaussianNoiseModel
        c = Curve([1., 2., 3.], [4., 5., 6.])
        nm = MvGaussianNoiseModel(1., 3., 5.)
        noise = add_noise(c, nm)
        @test length(noise) == 3
        @test noise.y |> first |> x -> x.particles |> length == 10_000
        noise = add_noise(c, nm, 1000)
        @test noise.y |> first |> x -> x.particles |> length == 1000
    end
    let # use GaussianNoiseModel
        c = Curve([1., 2., 3.], [4., 5., 6.])
        nm = GaussianNoiseModel(1.)
        noise = add_noise(c, nm)
        @test length(noise) == 3
        @test noise.y |> first |> x -> x.particles |> length == 10_000
        noise = add_noise(c, nm, 1000)
        @test noise.y |> first |> x -> x.particles |> length == 1000
    end
end


@testset "generate_noise() (correlated)" begin
    let
        nm = MvGaussianNoiseModel(1.0, 3.0, 5.0)
        noise = generate_noise(nm, 100, 2000)
        @test length(noise) == 100
        @test noise |> first |> x -> x.particles |> length == 2000
        ns = NoiseSample(mci.get_draw.(1, noise))
        noise2 = generate_noise(ns, 100)
        @test noise2 |> first |> x -> x.particles |> length == 100
    end
end


@testset "generate_noise() (uncorrelated)" begin
    let
        nm = GaussianNoiseModel(1.0)
        noise = generate_noise(nm, 100, 2000)
        @test length(noise) == 100
        @test noise |> first |> x -> x.particles |> length == 2000
    end
end


@testset "fit_noise()" begin
    # this also tests estimate_autocov()
    let
        seed!(2)
        samples = 100
        len = 2_000
        for (δxᵢ, αᵢ, λᵢ) ∈ [(1.0, 3.0, 5.0),
                             (0.5, 2.0, 0.5),
                             (2.0, 0.5, 3.0)]
            nm = MvGaussianNoiseModel(δxᵢ, αᵢ, λᵢ)
            noise = generate_noise(nm, len, samples)
            A = Float64[]
            Λ = Float64[]
            for i in 1:samples
                x = δxᵢ:δxᵢ:(len * δxᵢ) |> collect
                ns = NoiseSample(x, mci.get_draw.(i, noise))
                noise_param = fit_noise(ns)
                push!(A, noise_param.α)
                push!(Λ, noise_param.λ)
            end
            # To approach the true values, we need to have very long
            # noise samples, which gets computationally expensive
            # quickly. As a workaround, shorter noise samples
            # are used here, and the fitted parameters are averaged
            # over several trials.
            # When using a noise sample length of 2000 and averaging
            # over 100 trials, the fitted parameters should reproduce
            # the input parameters with 5% accuray:
            @test mean(A) ≈ αᵢ atol=5e-2
            @test mean(Λ) ≈ λᵢ atol=5e-2
        end
    end
end