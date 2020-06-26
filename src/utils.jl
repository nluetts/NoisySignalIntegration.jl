"""
    crop(s::Curve{T}, left::T, right::T) where {T}

Crop a slice from a spectrum and return a new spectrum.
"""
function crop(s::Curve{T}, left::T, right::T) where {T}
    i = searchsortedlast(s.x, left)
    j = searchsortedfirst(s.x, right)
    return Curve(s.x[i:j], s.y[i:j])
end