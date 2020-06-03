"""
This file defines the type `Spectrum` and basic functions
to manipulate and plot a spectrum.
"""

"""
    Spectrum(x::Array{T<:Real, 1}, y::Array{T<:Real, 1})
"""
struct Spectrum{T<:Real}
    x::Array{T, 1}
    y::Array{T, 1}
end

# basic methods

Base.:(==)(s0::Spectrum, s1::Spectrum) = s0.x == s1.x && s0.y == s1.y

Base.length(s::Spectrum) = length(s.x)
Base.lastindex(s::Spectrum) = length(s)
Base.getindex(s::Spectrum, i::Integer) = (s.x[i], s.y[i])
Base.getindex(s::Spectrum, r::UnitRange) =Spectrum(s.x[r], s.y[r])
function Base.iterate(s::Spectrum, i::Int64)
    if i > length(s)
        return nothing
    else
        return s[i], i + 1
    end
end
Base.iterate(s::Spectrum) = iterate(s, 1)
Base.IteratorSize(itertype::Type{Spectrum}) = Base.HasLength()

Base.vcat(s0::Spectrum, s1::Spectrum) = Spectrum(vcat(s0.x, s1.x), vcat(s0.y, s1.y))

function Base.show(io::IO, spectrum::Spectrum{T}) where {T}
    println("Spectrum{$T}, $(length(spectrum)) datapoints")
    if length(spectrum) > 10
        ps = [spectrum[1:5]; spectrum[end-4:end]]
        indicator = "    ⋮\n"
    else
        ps = spectrum
        indicator = ""
    end
    for (i, p) in enumerate(ps)
        @printf "    %-5e %-5e\n" p...
        i == 5 ? print(indicator) : nothing
    end
end

Plots.plot!(p::Plots.Plot, s::Spectrum, args...; kw...) = plot!(p, s.x, s.y, args...; kw...)
Plots.plot(s::Spectrum, args...; kw...) = plot!(plot(), s.x, s.y, args...; kw...)

"""
    function crop(slc::Slice) ::Spectrum

Crop a slice from a spectrum and return a new spectrum.
"""
function crop(s::Spectrum{T}, left::T, right::T) where {T}
    i = searchsortedlast(s.x, left)
    j = searchsortedfirst(s.x, right)
    return Spectrum(s.x[i:j], s.y[i:j])
end