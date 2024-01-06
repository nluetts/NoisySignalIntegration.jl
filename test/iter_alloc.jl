module IterAllocTest
using NoisySignalIntegration

@inline function lininterp(x, x₀, x₁, y₀, y₁)
    δx = x₁ - x₀
    y = (y₁ * (x - x₀) + y₀ * (x₁ - x)) / δx
    return y
end

const CurveIterState = Tuple{Int64,Bool}

struct MyVec
    x::Vector{Float64}
    y::Vector{Float64}
end

function Base.iterate(
    iter::Tuple{MyVec,Float64,Float64},
    state::Tuple{Int64,Bool}
)::Union{Nothing,Tuple{Tuple{Float64,Float64},Tuple{Int64,Bool}}}
    mv, left, right = iter
    index, return_right = state

    while index < length(mv.x)
        xi = mv.x[index]
        xj = mv.x[index+1]
        yi = mv.y[index]
        yj = mv.y[index+1]

        if xj <= left
            # not yet in bounds, continue while loop
            index += 1
            continue
        elseif xi >= right
            # moved outside bounds, stop iteration
            return nothing
        elseif xi >= left && xj <= right
            # xi and xj fall within bounds, return data point
            return ((xi, yi), (index + 1, false))
        elseif xi < left < xj && right >= xj
            # left bound falls between xi and xj,
            # interpolate and return data point at `left`
            yleft = lininterp(left, xi, xj, yi, yj)
            return ((left, yleft), (index + 1, false))
        elseif xi < right < xj && left <= xi
            # right bound falls between xi and xj
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
        elseif xi < left <= right < xj
            # very rare case where `left` and `right` fall
            # both between xi and xj
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
    return nothing
end

function Base.iterate(
    iter::Tuple{MyVec,Float64,Float64},
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


function main()
    x = range(start=-1, stop=10, length=1000) |> collect
    y = zeros(Float64, 1000)
    mv = MyVec(x, y)
    @time test_iter(mv)
    @time test_trapz(mv)
end

end
