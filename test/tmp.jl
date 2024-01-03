module New

import Base.iterate

using Distributions
using NoisySignalIntegration
using NoisySignalIntegration: AbstractNoiseModel, lininterp
using Plots

mutable struct LeftRightIterator{T}
    curve::Curve{T}
    left::T
    right::T
    index::Integer
    next::Union{Nothing,Tuple{T,T}}
    depleted::Bool
end

function iterate_left_right(crv::Curve{T}, left::T, right::T) where {T}
    left, right = min(left, right), max(left, right)
    depleted = false
    if crv.x[end] <= left || crv.x[1] >= right || left == right
        depleted = true
    end
    LeftRightIterator{T}(crv, left, right, 1, nothing, depleted)
end

function Base.iterate(lri::LeftRightIterator{T})::Union{Nothing,Tuple{Tuple{T,T},LeftRightIterator}} where {T}
    lri.depleted && return nothing
    if !isnothing(lri.next)
        lri.depleted = true
        return (lri.next, lri)
    end

    while true
        x1 = get(lri.curve.x, lri.index, nothing)
        x2 = get(lri.curve.x, lri.index + 1, nothing)
        y1 = get(lri.curve.y, lri.index, nothing)
        y2 = get(lri.curve.y, lri.index + 1, nothing)

        lri.index += 1

        # get state
        state = let
            left_between_x1x2 = (isnothing(x1) || isnothing(x2)) ? false : x1 < lri.left < x2
            right_between_x1x2 = (isnothing(x1) || isnothing(x2)) ? false : x1 < lri.right < x2
            in_bounds_x1 = isnothing(x1) ? nothing : lri.left <= x1 <= lri.right
            in_bounds_x2 = isnothing(x2) ? nothing : lri.left <= x2 <= lri.right
            (in_bounds_x1, in_bounds_x2, left_between_x1x2, right_between_x1x2)
        end

        # common case where we are outside the bounds
        # we have two sub-cases to deal with:
        # x2 < a -> continue (is executed implicitly by doing nothing)
        # x1 > b -> depleted
        if state == (false, false, false, false)
            if x1 > lri.right
                @show lri.depleted = true
                return nothing
            end
            # common case where we are inside the bounds
            # but did not just enter them or are about to exit them
        elseif state == (true, true, false, false)
            return ((x1, y1), lri)
            # moving inside the bounds
            # x1   a   x2   b -> return (a, y at a)
        elseif state == (false, true, true, false)
            return ((lri.left, lininterp(lri.left, x1, x2, y1, y2)), lri)
            # exiting the bounds
            # a   x1   b   x2 -> return (x1, y1)
            # -> set next_b = Some((b, y at b))
        elseif state == (true, false, false, true)
            lri.next = (lri.right, lininterp(lri.right, x1, x2, y1, y2))
            return ((x1, y1), lri)
            # no elements left -> depleted
        elseif state == (nothing, nothing, false, false)
            lri.depleted = true
            return nothing
            # only one element left, but not inside bounds
            # -> depleted (4 states)
        elseif !state[1] && isnothing(state[2])
            lri.depleted = true
            return nothing
            # only one element left, which is inside bounds
            # -> return (x1, y1) -> depleted (4 states)
        elseif state[1] && isnothing(state[2])
            lri.depleted = true
            return ((x1, y1), lri)
            # rare case where x2 falls on a
            # x1   x2 = a   b -> continue
        elseif state == (false, true, false, false)
            continue
            # rare case where x1 falls on b
            # a   x1 = b   x2 -> return (x1, y1)
        elseif state == (true, false, false, false)
            return ((x1, y1), lri)
            # very rare case where x1   a   b   x2
            # -> return (a, y at a)
            # -> set next_b = Some((b, y at b))
            # (note that a == b is not possible,
            # because this case is explicitly handled
            # when the iterator is constructed)
        elseif state == (false, false, true, true)
            self.next = (lri.right, lininterp(lri.right, x1, x2, y1, y2))
            return ((lri.left, lininterp(lri.left, x1, x2, y1, y2)), lri)
            # unreachable states
        else
            throw("Reached an iterator state that should be impossible, consider writing a bug report!")
        end
    end
end

Base.iterate(lri::LeftRightIterator, _::LeftRightIterator) = Base.iterate(lri::LeftRightIterator)

function update!(lri::LeftRightIterator{T}, crv::Curve{T}) where {T}
    length(crv) != length(lri.curve) && throw("Can only update iterator with curve of same length!")
    for i in 1:length(lri.curve)
        lri.curve.x[i] = crv.x[i]
        lri.curve.y[i] = crv.y[i]
    end
    lri.next = nothing
    lri.index = 1
    lri.depleted = false
end

function main()
    spectrum = NoisySignalIntegration.testdata_1()
    slice_bands = crop(spectrum, 5.0, 40.0)
    slice_noise = crop(spectrum, 40.0, 100.0)
    noise_model = fit_noise(NoiseSample(slice_noise, 3))
    ucrv = add_noise(slice_bands, noise_model)
    cov = get_cov(noise_model, length(slice_bands))
    iter = iterate_left_right(slice_bands, 10.37, 20.52)
    for i in 1:10000
        crv = NoisySignalIntegration.get_draw(i, ucrv)
        update!(iter, crv)
        for (x1, y1) in iter
            # print(x1, " ", y1, "\r")
        end
    end
end

end
