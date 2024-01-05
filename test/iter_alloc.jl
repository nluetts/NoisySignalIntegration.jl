module IterAllocTest
using NoisySignalIntegration

@inline function lininterp(x, x₀, x₁, y₀, y₁)
    δx = x₁ - x₀
    y = (y₁ * (x - x₀) + y₀ * (x₁ - x)) / δx
    return y
end

struct MyVec
    x::Vector{Float64}
    y::Vector{Float64}
end

function Base.iterate(
    iter::Tuple{MyVec,Float64,Float64},
    state::Tuple{UnitRange{Int64},UnitRange{Int64},Int64}
)::Union{Nothing,Tuple{Tuple{Float64,Float64},Tuple{UnitRange{Int64},UnitRange{Int64},Int64}}}
    mv, left, right = iter
    left_position, right_position, index = state

    if index > right_position.start || index > length(mv.x)
        return nothing
    end
    if index >= left_position.start && index < right_position.stop
        return ((mv.x[index], mv.y[index]), (left_position, right_position, index + 1))
    elseif index < left_position.stop
        yl = lininterp(left, mv.x[left_position.stop], mv.x[left_position.start], mv.y[left_position.stop], mv.y[left_position.start])
        return ((left, yl), (left_position, right_position, left_position.start))
    elseif index == right_position.start
        yr = lininterp(right, mv.x[right_position.stop], mv.x[right_position.start], mv.y[right_position.stop], mv.y[right_position.start])
        return ((right, yr), (left_position, right_position, index + 1))
    else
        return nothing
    end
end

function Base.iterate(
    iter::Tuple{MyVec,Float64,Float64},
)::Union{Nothing,Tuple{Tuple{Float64,Float64},Tuple{UnitRange{Int64},UnitRange{Int64},Int64}}}

    mv, left, right = iter
    left, right = min(left, right), max(left, right)
    left_position = searchsorted(mv.x, left)
    right_position = searchsorted(mv.x, right)
    Base.iterate(iter, (left_position, right_position, 1))
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
