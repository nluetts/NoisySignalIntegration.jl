"""
This file defines the type `Spectrum` and basic functions
to manipulate and plot a spectrum.
"""

"""
    AbstractSpectrum{T<:AbstractFloat}

Abstract supertype for all spectral data. Subtypes need to implement fields `x` and `y`
(xy-data of spectrum).
"""
abstract type AbstractSpectrum{T<:AbstractFloat} end

Base.length(s::AbstractSpectrum) = length(s.x)
Base.lastindex(s::AbstractSpectrum) = length(s)
Base.getindex(s::AbstractSpectrum, i::Integer) = (s.x[i], s.y[i])
Base.getindex(s::AbstractSpectrum, r::UnitRange) =Spectrum(s.x[r], s.y[r])
Base.iterate(s::AbstractSpectrum, i::Int64) = i > length(s) ? nothing : (s[i], i + 1)
Base.iterate(s::AbstractSpectrum) = iterate(s, 1)
Base.IteratorSize(itertype::Type{AbstractSpectrum}) = Base.HasLength()

function Base.show(io::IO, spectrum::AbstractSpectrum{T}) where {T}
    println("$(typeof(spectrum)), $(length(spectrum)) datapoints")
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

Plots.plot!(p::Plots.Plot, s::AbstractSpectrum, args...; kw...) = plot!(p, s.x, s.y, args...; kw...)
Plots.plot(s::AbstractSpectrum, args...; kw...) = plot!(plot(), s.x, s.y, args...; kw...)

"""
    Spectrum(x::Array{T<:AbstractFloat, 1}, y::Array{T<:AbstractFloat, 1})
"""
struct Spectrum{T} <: AbstractSpectrum{T}
    x::Vector{T}
    y::Vector{T}
    function Spectrum(x::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
        length(x) == length(y) || throw(ArgumentError("x and y need to have the same length."))
        return new{T}(x, y)
    end
end

# basic methods

Base.:(==)(s0::Spectrum, s1::Spectrum) = s0.x == s1.x && s0.y == s1.y
Base.:(+)(s::Spectrum{T}, y::Vector{T}) where {T}= Spectrum(s.x, s.y .+ y)
Base.vcat(s0::Spectrum, s1::Spectrum) = Spectrum(vcat(s0.x, s1.x), vcat(s0.y, s1.y))

"""
    function crop(slc::Slice) ::Spectrum

Crop a slice from a spectrum and return a new spectrum.
"""
function crop(s::Spectrum{T}, left::T, right::T) where {T}
    i = searchsortedlast(s.x, left)
    j = searchsortedfirst(s.x, right)
    return Spectrum(s.x[i:j], s.y[i:j])
end