using Test
using MCIntegrate

@testset "Curve" begin

    curve = Curve([3, 4, 5])
    @test curve == Curve([1, 2, 3], [3, 4, 5])

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

    # x and y need to have same length
    @test_throws ArgumentError Curve([1.0], [4.0, 4.0])
    @test_throws ArgumentError Curve([1.0, 2.0], [4.0])
end;

@testset "create Curves from collections" begin
    x = 1:10
    y = 11:20
    c = Curve(x, y)
    @test eltype(c) == Int
    @test c.x == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    @test c.y == [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    c = Curve(1:3, [1. , 2., 3.])
    @test eltype(c) == Float64
    @test c.x == [1., 2., 3.]
    @test c.y == [1., 2., 3.]
    c = Curve([1. , 2., 3.], 1:3)
    @test eltype(c) == Float64
    @test c.x == [1., 2., 3.]
    @test c.y == [1., 2., 3.]
    # element types must be subtypes of Number
    @test_throws MethodError Curve(["a", "c"], ["b", "d"])
end


@testset "crop()" begin
    c = Curve(1:10, 11:20)
    @test crop(c, 4, 6) == Curve([4, 5, 6], [14, 15, 16])
    c = Curve(1:10, collect(Float64, 11:20))
    @test crop(c, 3.1, 6.9) == Curve([4., 5., 6.], [14., 15., 16.])
    @test crop(c, 3.99999999999999, 6.00000000000001) == Curve([4., 5., 6.], [14., 15., 16.])
    @test crop(c, 3.1, 6.9) == Curve([4., 5., 6.], [14., 15., 16.])
    @test crop(c, 3.1, 7.0) == Curve([4., 5., 6., 7.], [14., 15., 16., 17.])
    @test crop(c, 3.0, 7.0) == Curve([3., 4., 5., 6., 7.], [13., 14., 15., 16., 17.])
end