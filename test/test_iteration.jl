@testset "test left-right iterator" begin
    ds = let
        x = -1:10 |> collect
        Curve{Float64}(x, x * 2)
    end
    # iterate window inside data range
    # both left and right fall on data points
    iter = Iterators.Stateful((ds, 1.0, 3.0))
    @test iterate(iter) |> first == (1.0, 2.0)
    @test iterate(iter) |> first == (2.0, 4.0)
    @test iterate(iter) |> first == (3.0, 6.0)
    @test iterate(iter) |> isnothing
    # right falls on data point, left does not
    iter = Iterators.Stateful((ds, 1.5, 3.0))
    @test iterate(iter) |> first == (1.5, 3.0)
    @test iterate(iter) |> first == (2.0, 4.0)
    @test iterate(iter) |> first == (3.0, 6.0)
    @test iterate(iter) |> isnothing
    # left falls on data point, right does not
    iter = Iterators.Stateful((ds, 1.0, 2.5))
    @test iterate(iter) |> first == (1.0, 2.0)
    @test iterate(iter) |> first == (2.0, 4.0)
    @test iterate(iter) |> first == (2.5, 5.0)
    @test iterate(iter) |> isnothing
    # left and right do not fall on data point
    iter = Iterators.Stateful((ds, 1.5, 3.4))
    @test iterate(iter) |> first == (1.5, 3.0)
    @test iterate(iter) |> first == (2.0, 4.0)
    @test iterate(iter) |> first == (3.0, 6.0)
    @test iterate(iter) |> first == (3.4, 6.8)
    @test iterate(iter) |> isnothing
    # left to right interval is smaller than x-spacing,
    # left falls on x-grid
    iter = Iterators.Stateful((ds, 1.0, 1.5))
    @test iterate(iter) |> first == (1.0, 2.0)
    @test iterate(iter) |> first == (1.5, 3.0)
    @test iterate(iter) |> isnothing
    # left to right interval is smaller than x-spacing,
    # right falls on x-grid
    iter = Iterators.Stateful((ds, 1.5, 2.0))
    @test iterate(iter) |> first == (1.5, 3.0)
    @test iterate(iter) |> first == (2.0, 4.0)
    @test iterate(iter) |> isnothing
    # left to right interval is exactly as wide as x-spacing
    iter = Iterators.Stateful((ds, 1.0, 2.0))
    @test iterate(iter) |> first == (1.0, 2.0)
    @test iterate(iter) |> first == (2.0, 4.0)
    @test iterate(iter) |> isnothing
    # iterate starting outside data range
    iter = Iterators.Stateful((ds, -2.0, 2.0))
    @test iterate(iter) |> first == (-1.0, -2.0)
    @test iterate(iter) |> first == (0.0, 0.0)
    @test iterate(iter) |> first == (1.0, 2.0)
    @test iterate(iter) |> first == (2.0, 4.0)
    @test iterate(iter) |> isnothing
    # starting further outside bound should not matter
    iter = Iterators.Stateful((ds, -20.0, 2.0))
    @test iterate(iter) |> first == (-1.0, -2.0)
    # iterate ending outside data range
    iter = Iterators.Stateful((ds, 8.0, 11.0))
    @test iterate(iter) |> first == (8.0, 16.0)
    @test iterate(iter) |> first == (9.0, 18.0)
    @test iterate(iter) |> first == (10.0, 20.0)
    @test iterate(iter) |> isnothing
    # iterate starting and ending outside data range
    iter = Iterators.Stateful((ds, -8.0, 18.0))
    @test iterate(iter) |> first == (-1.0, -2.0)
    [iterate(iter) for _ in 1:10] # skip central elements
    @test iterate(iter) |> first == (10.0, 20.0)
    @test iterate(iter) |> isnothing
    # left and right fall between data points
    iter = Iterators.Stateful((ds, 2.25, 2.75))
    @test iterate(iter) |> first == (2.25, 4.5)
    @test iterate(iter) |> first == (2.75, 5.5)
    @test iterate(iter) |> isnothing
end
