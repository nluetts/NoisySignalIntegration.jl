"""
    crop(s::Curve, left, right)

Crop a slice from a curve and return a new curve.
"""
function crop(s::Curve, left::Float64, right::Float64)
    i = searchsortedlast(s.x, left)
    j = searchsortedfirst(s.x, right)
    return Curve(s.x[i:j], s.y[i:j])
end
crop(s::Curve, left, right) = crop(s, Float64(left), Float64(right))

"""
    clone(wb::WidthBound, loc::Float64)

Create a clone of a WidthBound that retrieves samples from the parent object.

This function returns a WidthBoundClone object.
A WidthBound together with a WidthBoundClone is used to integrate two
parts of the curve with the same width in each Monte-Carlo draw.
"""
function clone(wb::WidthBound, loc::Float64)
    return WidthBoundClone(Float64(loc), wb)
end