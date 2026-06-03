struct UncertainPseudoVoigtFit{T,N}
    area::Particles{T,N}
    center::Particles{T,N}
    width::Particles{T,N}
    mixing::Particles{T,N}
    offset::Particles{T,N}
    slope::Particles{T,N}
end

struct PseudoVoigtFit{T}
    area::T
    center::T
    width::T
    mixing::T
    offset::T
    slope::T
end

function get_draw(n, f::UncertainPseudoVoigtFit{T,N})::PseudoVoigtFit{T} where {T,N}
    PseudoVoigtFit(f.area, f.center, f.width, f.mixing, f.offset, f.slope)
end

function Base.show(io::IO, ::MIME"text/plain", obj::UncertainPseudoVoigtFit)
    # Print type header
    println(io, "UncertainPseudoVoigtFit{$(eltype(obj.area)),$(size(obj.area, 2))}")

    # Define fields and their names
    fields = (:area, :center, :width, :mixing, :offset, :slope)
    # Calculate max field name length for alignment
    max_len = maximum(length(string(f)) for f in fields)

    for f in fields
        val = getfield(obj, f)
        m = pmean(val)
        s = pstd(val)
        # Format: field_name = mean ± std
        # Use Printf for alignment and formatting
        @printf(io, " %s = %6.3g ± %6.2g\n", lpad(string(f), max_len), m, s)
    end
end

function pvoigt_profile(x, area, center, width, mixing, offset, slope)
    u = (x - center) / (width / 2)

    # The Lorentzian fraction is calculated via a sigmoid function
    # instead of capping the input mixing ratio to the range [0, 1].
    # This is done for numerical stability when fitting.
    # The factor 7 is chosen so that for the `mixing = -1` the mixing
    # ratio is 99.9% Gaussian.
    ratio = 1 / (1 + exp(-7 * mixing))

    gaussian_term = (1 - ratio) * √log(2) / (width * √pi) * exp(-log(2) * u^2)
    lorentzian_term = ratio / (pi * width * (1 + u^2))

    return area * (gaussian_term + lorentzian_term) + offset + x * slope
end

pvoigt_profile(x, params) = pvoigt_profile.(x, params...)

function pvoigt_profile(xs::S, pvoigt_fit::PseudoVoigtFit{T}) where {S<:AbstractVector,T}
    f = pvoigt_fit
    ys = zeros(length(xs))
    for (i, x) in enumerate(xs)
        ys[i] += pvoigt_profile(x, f.area, f.center, f.width, f.mixing, f.offset, f.slope)
    end
    ys
end

"""Height of Pseudo-Voigt peak."""
function pvoigt_peak(area, width, mixing_fraction)
    area * ((1 - mixing_fraction) * √log(2) / (width * √pi) + (mixing_fraction / (width * pi)))
end


function fit_pvoigt(
    curve::Curve{T},
    left::T,
    right::T,
) where {T<:AbstractFloat}
    mask = curve.x .> left .&& curve.x .<= right
    xs = curve.x[mask]
    ys = curve.y[mask]

    # Guess
    guess = let
        #! format: off
        mixing   = 0.0 # 50-50 Gaussian to Lorentzian
        height   = maximum(ys) - minimum(ys)
        center   = (right + left) * 0.5
        width    = abs(right - left) * 0.25 # 1/4 of the fit-window
        area     = height / pvoigt_peak(1.0, width, mixing)
        offset   = 0.0
        slope    = 0.0
        #! format: on

        [area, center, width, mixing, offset, slope]
    end

    params = curve_fit(pvoigt_profile, xs, ys, guess) |> coef
    PseudoVoigtFit(params...)
end

function mc_fit(
    uc::UncertainCurve{Float64,N},
    bnds::Vector{UncertainBound{Float64,M}};
)::Vector{UncertainPseudoVoigtFit{Float64,N}} where {M,N}

    # TODO: This could be a warning and we just take the smaller amount of samples
    M != N && error("Samples sizes of bounds and uncertain curve incompatible ($N != $M)")

    pfits = Array{PseudoVoigtFit{Float64}}(undef, N, length(bnds))
    for i ∈ 1:N
        i % 1000 == 0 && print("Processing draw $i/$N \r")
        cᵢ = get_draw(i, uc)
        for (j, b) in enumerate(bnds)
            xₗ, xᵣ = get_draw(i, b)
            pfits[i, j] = fit_pvoigt(cᵢ, xₗ, xᵣ)
        end
    end

    # Pack results into UncertainPseudoVoigtFit
    fit_results = []
    buf = Array{Float64}(undef, N, 6)
    for j in 1:length(bnds)
        for i in 1:N
        #! format: off
            fit = pfits[i, j]
            buf[i, 1] = fit.area
            buf[i, 2] = fit.center
            buf[i, 3] = fit.width
            buf[i, 4] = fit.mixing
            buf[i, 5] = fit.offset
            buf[i, 6] = fit.slope
        end
        area   = Particles(buf[:, 1])
        center = Particles(buf[:, 2])
        width  = Particles(buf[:, 3])
        mixing = Particles(buf[:, 4])
        offset = Particles(buf[:, 5])
        slope  = Particles(buf[:, 6])
        #! format: on
        push!(fit_results, UncertainPseudoVoigtFit(area, center, width, mixing, offset, slope))
    end
    fit_results
end
