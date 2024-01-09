module IterAllocTest
using NoisySignalIntegration
using Test

@inline function lininterp(x, x₀, x₁, y₀, y₁)
    δx = x₁ - x₀
    y = (y₁ * (x - x₀) + y₀ * (x₁ - x)) / δx
    return y
end

function eltype(::Tuple{Curve{Float64},Float64,Float64})
    return Tuple{Float64,Float64}
end

function Base.iterate(
    iter::Tuple{Curve{Float64},Float64,Float64},
    state::Tuple{Int64,Bool}
)::Union{Nothing,Tuple{Tuple{Float64,Float64},Tuple{Int64,Bool}}}
    mv, left, right = iter
    index, return_right = state

    while index < length(mv.x)
        xi = mv.x[index]
        xj = mv.x[index+1]
        yi = mv.y[index]
        yj = mv.y[index+1]

        state = let
            in_bounds_xi = left <= xi <= right
            in_bounds_xj = left <= xj <= right
            left_between_ij = xi < left < xj
            right_between_ij = xi < right < xj
            in_bounds_xi, in_bounds_xj, left_between_ij, right_between_ij
        end

        if state == (false, false, false, false)
            # common case where we are outside the bounds
            # we have two sub-cases to deal with:
            # xj < left -> continue
            # xi > right -> depleted
            if xi > right
                return nothing
            end
            # not yet in bounds, continue loop
            index += 1
            continue
        elseif state == (true, true, false, false)
            # common case where we are inside the bounds
            # but did not just enter them or are not just about to exit them
            return ((xi, yi), (index + 1, false))
        elseif state == (false, true, true, false)
            # moving inside the bounds
            # xi   left   xj   right -> return (left, y at left)
            return ((left, lininterp(left, xi, xj, yi, yj)), (index + 1, false))
        elseif state == (true, false, false, true)
            # exiting the bounds
            # left   xi   right   xj -> return (xi, yj),
            # then return (right, y at right) in next iteration
            if return_right
                # if we already returned xi for this index,
                # interpolate and return data point at `right` 
                yright = lininterp(right, xi, xj, yi, yj)
                return ((right, yright), (index + 1, false))
            else
                # we have to return (xi, yi) before returning
                # interpolated data point at `right`;
                # note that index is _not_ incremented and
                # `return_right` flag is set
                return ((xi, yi), (index, true))
            end
        elseif state == (false, true, false, false)
            # rare case where xj falls on left
            # xi   xj = left   right -> continue
            index += 1
            continue
        elseif state == (true, false, false, false)
            # rare case where xi falls on right
            # left   xi = right   xj -> return (xi, yi)
            return ((xi, yi), (index + 1, false))
        elseif state == (false, false, true, true)
            # very rare case where xi   left   right   xj
            # -> return (left, y at left) and (right, y at right) in next iteration
            # (note that left == right is not possible,
            # because this case is explicitly handled
            # when the iterator is constructed)
            if return_right
                yright = lininterp(right, xi, xj, yi, yj)
                return ((right, yright), (index + 1, false))
            else
                yleft = lininterp(left, xi, xj, yi, yj)
                return ((left, yleft), (index, true))
            end
        else
            throw("Met an iterator state that should not have been possible")
        end
    end
    # check if last element is still within bounds
    if index == length(mv.x) && left <= mv.x[index] <= right
        return ((mv.x[index], mv.y[index]), (index + 1, false))
    else
        return nothing
    end
end

function Base.iterate(
    iter::Tuple{Curve{Float64},Float64,Float64},
)::Union{Nothing,Tuple{Tuple{Float64,Float64},Tuple{Int64,Bool}}}

    mv, left, right = iter
    left, right = min(left, right), max(left, right)
    if left == right || mv.x[end] <= left || mv.x[1] >= right
        return nothing
    end
    Base.iterate(iter, (1, false))
end

function test_iter(mv)
    A = 0.0
    for _ in 1:10000
        for (xi, yi) in (mv, 5.5, 7.8)
            A += xi + yi
        end
    end
    A
end

function test_trapz(mv)
    A = 0.0
    for _ in 1:10000
        A += NoisySignalIntegration.trapz(mv.x, mv.y, 5.5, 7.8)
    end
    A
end

function default(::Type{Curve{Float64}})
    x = -1:10 |> collect
    Curve{Float64}(x, x * 2)
end


function main()
    mv = default(Curve{Float64})
    @time test_iter(mv)
    @time test_trapz(mv)
end

function run_tests()
    @testset "test left-right iterator" begin
        ds = default(Curve{Float64})
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
end

end
