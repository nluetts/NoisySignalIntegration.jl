using Test
using MCIntegrate

@testset "AbstractCurve" begin

    curve = Curve([1.0, 2.0], [3.0, 4.0])

    # test basic operations
    @test length(curve) == 2
    @test curve[1] == (1.0, 3.0)
    @test curve[end] == (2.0, 4.0)
    @test [p for p in curve] == [(1.0, 3.0), (2.0, 4.0)] # test iteration
    @test [curve; curve] == Curve([1.0, 2.0, 1.0, 2.0], [3.0, 4.0, 3.0, 4.0])
    @test curve + curve.y == Curve([1.0, 2.0], [6.0, 8.0])
    @test minimum(curve) == (1.0, 3.0)
    @test maximum(curve) == (2.0, 4.0)
    @test curve + 2.0 == Curve([1.0, 2.0], [5.0, 6.0])
    @test curve - 2.0 == Curve([1.0, 2.0], [1.0, 2.0])
    @test curve * 2.0 == Curve([1.0, 2.0], [6.0, 8.0])
    @test curve / 2.0 == Curve([1.0, 2.0], [1.5, 2.0])
    @test 2.0 + curve == Curve([1.0, 2.0], [5.0, 6.0])
    @test 2.0 - curve == Curve([1.0, 2.0], [-1.0, -2.0])
    @test 2.0 * curve == Curve([1.0, 2.0], [6.0, 8.0])
    @test 2.0 / curve == Curve([1.0, 2.0], [2.0/3.0, 0.5])

    noise = Noise([1.0, 2.0], [4.0, 4.0])
    @test length(noise) == 2
    @test noise[end][1] == 2.0
    @test noise[end][2] + 1 ≈ 1.0
    @test [noise; noise] == Noise([1.0, 2.0, 1.0, 2.0], [4.0, 4.0, 4.0, 4.0])

    # x and y need to have same length
    @test_throws ArgumentError Curve([1.0], [4.0, 4.0])
    @test_throws ArgumentError Noise([1.0], [4.0, 4.0])
    @test_throws ArgumentError Curve([1.0, 2.0], [4.0])
    @test_throws ArgumentError Noise([1.0, 2.0], [4.0])
end;