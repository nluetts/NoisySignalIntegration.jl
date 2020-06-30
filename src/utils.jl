"""
    crop(s::Curve, left, right)

Crop a slice from a spectrum and return a new spectrum.
"""
function crop(s::Curve, left::Float64, right::Float64)
    i = searchsortedlast(s.x, left)
    j = searchsortedfirst(s.x, right)
    return Curve(s.x[i:j], s.y[i:j])
end
crop(s::Curve, left, right) = crop(s, Float64(left), Float64(right))

clone(wb::WidthBound, loc::Float64) = WidthBoundClone(Float64(loc), wb)