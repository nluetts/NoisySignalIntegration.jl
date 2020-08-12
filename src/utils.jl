"""
    crop(s::Curve, left, right)

Crop a slice from a curve and return a new curve.
"""
function crop(s::AbstractCurve, left::Float64, right::Float64)
    i = searchsortedlast(s.x, left)
    j = searchsortedfirst(s.x, right)
    return Curve(s.x[i:j], s.y[i:j])
end
crop(s::AbstractCurve, left, right) = crop(s, Float64(left), Float64(right))