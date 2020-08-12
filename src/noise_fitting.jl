"""
    gauss_kernel(Δx, β)

Fit function for noise autocorrelation.

    β = [α, λ]
    f(δx) = α^2 * exp(-0.5 * (Δx/λ)^2)
"""
function gauss_kernel(Δx, β)
    α, λ = β
    return @. α^2 * exp(-0.5 * (Δx/λ)^2)
end

"""
    get_cov(δx::Float64, n::Integer, α::Float64, λ::Float64)

Return covariance matrix `Σ` with size `n` × `n` for the exponentiated quadratic kernel
with amplitude `α` and lag `λ` calculated on a grid of spacing `δx`.
"""
function get_cov(δx::Float64, n::Integer, α::Float64, λ::Float64)
    β = [α, λ]
    Σ = Array{Float64}(undef, n, n)
    for i = (1:n), j = (1:n)
        Δx = δx * (i - j)
        Σ[i,j] = gauss_kernel(Δx, β)
        if i == j
            Σ[i,j] *= 1.0000001 # for numerical stability (?)
        end
    end
    return Σ
end


function estimate_autocov(n::Curve)
    # check that spectral grid is equally spaced
    if !allapproxequal(diff(n.x))
        throw(ArgumentError("Noise analysis only works on an evenly spaced spectral grid."))
    end
    δx = abs(n.x[1]-n.x[2]) # x-values must be evenly spaced
    lag_units = begin
        N = convert(Integer, round(length(n)/2))
        collect(1:N)
    end
    lags = δx .* lag_units
    acov = autocov(n.y, lag_units)
    return lags, acov
end

"""
    function fit_noise(
        lags::Vector{Float64},
        autocov::Vector{Float64};
        α_guess=1.0,
        λ_guess=1.0
    )

Estimates parameters α and λ of the exponentiated quadratic covariance function
k(Δx) = α² exp(-0.5 (Δx/λ)²) fitted to autocovariance (k) vs. lag (Δx).
If a `Noise` instance is passes as first argument, its autocovarince is
estimated and passed to the fitting function.

Note: Change initial guess if fit does not converge to sensible result.
"""
function fit_noise(
    lags::Vector{Float64},
    autocov::Vector{Float64};
    α_guess=1.0,
    λ_guess=1.0
)
    fit = curve_fit(gauss_kernel, lags, autocov, [α_guess, λ_guess])
    return fit.param
end

function fit_noise(n::Curve; α_guess=1.0, λ_guess=1.0, detrend_order=2)
    lags, acov = detrend(n, detrend_order) |> estimate_autocov
    fit_noise(lags, acov; α_guess=α_guess, λ_guess=λ_guess)
end

"""
Generate correlated noise from a noise sample.
"""
function correlated_noise(noise_sample::Curve, len::Int, samples::Int; kw...)
    n = noise_sample
    α, λ = fit_noise(n; kw...)
    δx = n.x[2] - n.x[1]
    Σ = get_cov(δx, len, α, λ)
    return Particles(samples, MvNormal(Σ))
end
