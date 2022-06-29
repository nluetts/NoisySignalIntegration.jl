using MonteCarloMeasurements: pstd


num_part(::UncertainBound{T, N}) where {T, N} = N
ub_eltype(::UncertainBound{T, N}) where {T, N} = T


@testset "UncertainBound (default constructor)" begin
    ub = UncertainBound(1..2, 3..4)
    @test num_part(ub) == NP
    @test ub.left.particles |> length == NP
end


@testset "get_draw()" begin
    ps = Particles([1, 2, 3])
    ub = UncertainBound(ps, ps)
    for i ∈ 1:3
        @test nsi.get_draw(i, ub) == [i, i]
    end
end


@testset "UncertainBound (construct from distributions)" begin
    dleft = Normal(0, 1)
    dright = Uniform(2, 3)
    seed!(1)
    ub = UncertainBound(dleft, dright)
    @test num_part(ub) == 10_000 # default no. of samples from constructor
    @test ub_eltype(ub) == Float64  # Distributions.jl should yield Float64 ny default
    @test ub.left |> pstd ≈ 1 atol = 1e-4
    @test ub.right |> pstd ≈ 1/√(12) atol = 1e-4
    ub = UncertainBound(dleft, dright, 1000)
    @test num_part(ub) == 1000
end


@testset "UncertainBound (construct from position(s) and width)" begin
    
    seed!(1)

    # build an uncertain curve with known position of peaks
    # in each draw
    uc = begin
        x = collect(Float64, 0:0.1:20)
        y = zeros(Float64, 9, length(x))
        # make peaks fall in x-range 10:10.9 and 15:15.9
        # draw 1: peak at 10.0 and 15.0
        # draw 2: peak at 10.1 and 15.1
        # draw 3: peak at 10.2 and 15.2
        # etc.
        for i ∈ 1:9
            y[i, i + 100] = 1.
            y[i, i + 150] = 1.
        end
        UncertainCurve(x, Particles(y))
    end

    dwidth = Uniform(1, 1 + 1e-8) # I want essentially the same width everywhere for testing
    ub = UncertainBound(10.5, dwidth, uc)
    # with the above setup, the values of the `left` and `right` particles should be as follows
    # the deviation comes from the fact that the width distribution (Uniform) hast to have some small width
    @test ub.left.particles[1] ≈ 9.5
    @test ub.left.particles[5] ≈ 9.9
    @test ub.left.particles[9] ≈ 10.3
    @test ub.right.particles[1] ≈ 10.5
    @test ub.right.particles[5] ≈ 10.9
    @test ub.right.particles[9] ≈ 11.3

    # test creation of bounds with correlated width
    dwidth = Uniform(1, 2) # now the width should vary more, to test the correlation of width between the two generated bounds
    ub1, ub2 = UncertainBound([10.5, 15.5], dwidth, uc)
    for i ∈ 1:9
        @test ub1.left.particles[i] - ub1.right.particles[i] ≈ ub2.left.particles[i] - ub2.right.particles[i] atol = 1e-9
    end
end


@testset "UncertainBound creation throws error when distribution has invalid support region" begin
    seed!(1)
    uc = begin # create uncertain curve with one symmetric peak
        x = 0:0.1:10;
        y = @. exp(-(x - 5)^2)
        add_noise(Curve(x, y), GaussianNoiseModel(0.03))
    end;
    @test_throws ArgumentError UncertainBound(5.0, Normal(0, 1), uc)
end


@testset "scale_shift_beta()" begin
    seed!(1)
    β = scale_shift_beta(2, 2, 1, 2)
    xs = rand(β, 1_000_000)
    @test minimum(xs) ≈ 1 atol = 1e-3
    @test maximum(xs) ≈ 2 atol = 1e-3
    @test mean(xs) ≈ 1.5 atol = 1e-3
end
