FIT_MAX_TRIES = 10
FIT_SHUFFLE_GUESS = 0.05

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
    PseudoVoigtFit(
        f.area.particles[n],
        f.center.particles[n],
        f.width.particles[n],
        f.mixing.particles[n],
        f.offset.particles[n],
        f.slope.particles[n]
    )
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

    gaussian_term = (1 - mixing) * √log(2) / (width * √pi) * exp(-log(2) * u^2)
    lorentzian_term = mixing / (pi * width * (1 + u^2))

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
function pvoigt_peak(area, width, mixing)
    area * ((1 - mixing) * √log(2) / (width * √pi) + (mixing / (width * pi)))
end


function fit_pvoigt(
    curve::Curve{T},
    left::T,
    right::T;
    guess=nothing
) where {T<:AbstractFloat}
    mask = curve.x .> left .&& curve.x .<= right
    xs = curve.x[mask]
    ys = curve.y[mask]

    # Normalize data, simpler guess and (hopefully) better numeric stability
    xmin, xmean, xmax = minimum(xs), mean(xs), maximum(xs)
    xspan = abs(xmax - xmin)
    xt = @. (xs - xmean) / xspan
    ymin, ymean, ymax = minimum(ys), mean(ys), maximum(ys)
    yspan = abs(ymax - ymin)
    yt = @. (ys - ymin) / yspan

    # Note on normalizing the guess: Everything is rather simple, but
    # the slope and offset need some thinking:
    # 
    # If x', y', m' and b' are coordinates and parameters in the
    # transformed normalized coordinate system, we can transform via:
    # 
    # x' = (x - xmean) / xspan
    # y' = (y - ymin) / yspan
    # y' = m'x' + b'
    # (y - ymin) / yspan = m'(x - xmean) / xspan + b'
    # y = ymin + m' yspan/xspan (x - xmean) + b'yspan
    # → m = m'yspan/xspan
    # → b = ymin + b'yspan - m'yspan/xspan xmean
    # → b = ymin + b'yspan - m xmean
    # We can use these to transform the fitted parameters
    # or invert the equations to a normalize the guess:
    # m' = m xspan/yspan
    # b' = (b - ymin + m xmean)/yspan

    # Guess
    guess = isnothing(guess) ? let
        #! format: off
        mixing   = 0.5 # 50-50 Gaussian to Lorentzian
        height   = 1.0
        center   = 0.0
        width    = 0.2 # 1/3 of the fit-window
        area     = height / pvoigt_peak(1.0, width, mixing)
        offset   = 0.0
        # use start and end point as heuristic for baseline
        slope    = yt[end] - yt[1]
        #! format: on

        [area, center, width, mixing, offset, slope]
    end : let
        area, center, width, mixing, offset, slope = guess
        # Normalize guess
        area /= xspan * yspan
        center = (center - xmean) / xspan
        width /= xspan
        # mind that the slope was not updated yet!
        offset = (offset - ymin + slope * xmean) / yspan
        slope *= xspan / yspan

        [area, center, width, mixing, offset, slope]
    end

    # let
    #     plot(curve)
    #     area, center, width, mixing, offset, slope = guess
    #     area *= xspan * yspan
    #     center = center * xspan + xmean
    #     width *= xspan
    #     slope *= yspan / xspan
    #     # mind that slope was already updated!
    #     offset = ymin + yspan * offset - slope * xmean
    #     plot!(PseudoVoigtFit(area, center, width, mixing, offset, slope)) |> display
    # end

    # bounds
    lower = [
        0.0,    # area
        -0.5,   # center
        0.0,    # width
        0,      # mixing
        -Inf,   # offset
        -Inf,   # slope
    ]
    upper = [
        Inf,   # area
        0.5,   # center
        1.0,   # width
        1,     # mixing
        Inf,   # offset
        Inf,   # slope
    ]

    for k in 1:FIT_MAX_TRIES
        try
            res = curve_fit(pvoigt_profile, xt, yt, guess; upper=upper, lower=lower, show_trace=true)
            params = coef(res)
            area, center, width, mixing, offset, slope = params
            # Undo normalization
            area *= xspan * yspan
            center = center * xspan + xmean
            width *= xspan
            slope *= yspan / xspan
            # mind that slope was already updated!
            offset = ymin + yspan * offset - slope * xmean
            return PseudoVoigtFit(area, center, width, mixing, offset, slope)
        catch e
            @warn "Fitting failed with error: $(e)"
            @info "Current guess: $guess"
            @info "Retrying fit with randomized guess ... ($k/$FIT_MAX_TRIES)"
            # Randomize guess by adding a fraction of a standard deviation times the actual value
            guess = [g + randn() * FIT_SHUFFLE_GUESS * g for g in guess]
            if k == FIT_MAX_TRIES
                throw(error("Could not fit curve in 10 trials, bailing out."))
            end
        end
    end
end

function mc_fit(
    uc::UncertainCurve{Float64,N},
    bnds::Vector{UncertainBound{Float64,M}};
)::Vector{UncertainPseudoVoigtFit{Float64,N}} where {M,N}

    # TODO: This could be a warning and we just take the smaller amount of samples
    M != N && error("Samples sizes of bounds and uncertain curve incompatible ($N != $M)")

    pfits = Array{PseudoVoigtFit{Float64}}(undef, N, length(bnds))
    # We have to go bounds first here, because otherwise the guess that
    # is re-used from draw to draw will be applied to the wrong bound
    for (j, b) in enumerate(bnds)
        println("\nProcessing bound $j")
        guess = nothing
        for i ∈ 1:N
            i % 1000 == 0 && print("Processing draw $i/$N \r")
            cᵢ = get_draw(i, uc)
            xₗ, xᵣ = get_draw(i, b)
            f = fit_pvoigt(cᵢ, xₗ, xᵣ, guess=guess)
            guess = [f.area, f.center, f.width, f.mixing, f.offset, f.slope]
            pfits[i, j] = f
            if i == 10
                throw("Testing ...")
            end
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

mc_fit(uc::UncertainCurve{Float64,N}, bd::UncertainBound{Float64,M}) where {M,N} = mc_fit(uc, [bd])
