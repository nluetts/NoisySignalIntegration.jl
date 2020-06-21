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
    get_cov(δx::T, n::Integer, α::T, λ::T) where {T<:AbstractFloat}

Return covariance matrix `Σ` with size `n` × `n` for the exponentiated quadratic kernel
with amplitude `α` and lag `λ` calculated on a grid of spacing `δx`.
"""
function get_cov(δx::T, n::Integer, α::T, λ::T) where {T<:AbstractFloat}
    β = [α, λ]
    Σ = Array{T}(undef, n, n)
    for i = (1:n), j = (1:n)
        Δx = δx * (i - j)
        Σ[i,j] = gauss_kernel(Δx, β)
        if i == j
            Σ[i,j] *= 1.0000001 # for numerical stability (?)
        end
    end
    return Σ
end

#######################################################
### Sampling ##########################################
#######################################################

function sample(nm::GaussianNoiseModel{T}, m::Integer=100, n::Integer=1) where {T}
    return nm.σ * randn(T, m, n)
end


function sample(
    nm::MvGaussianNoiseModel{T},
    m::Integer=100,
    n::Integer=1
) where {T<:AbstractFloat}
    samples = rand(MvNormal(get_cov(nm.δx, m, nm.α, nm.λ)), n)
    n == 1 && return samples[:,1]
    return samples
end


function sample(b::LeftRightBound, samples::Integer=1)
    left = rand(b.left, samples)
    right = rand(b.right, samples)
    return hcat(left, right)
end


sample(wb::WidthBound, samples::Integer=1) = rand(wb.width, samples)
function sample(wb::WidthBound, s::Curve, samples::Integer=1) 
    ws = rand(wb.width, samples)
    spls = [left_right_from_peak(s.x, s.y, wb.loc, w) for w in ws]
    samples == 1 && return spls[1]
    return spls
end

#######################################################
### Fitting ###########################################
#######################################################

function estimate_autocov(n::Noise{T}) where {T<:AbstractFloat}
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
    fit_noise(lags::T, autocov::T; α_guess=1.0, λ_guess=1.0) where {T<:AbstractFloat}
    fit_noise(n::Noise{T}; α_guess=1.0, λ_guess=1.0) where {T}

Estimates parameters α and λ of the exponentiated quadratic covariance function
k(Δx) = α² exp(-0.5 (Δx/λ)²) fitted to autocovariance (k) vs. lag (Δx).
If a `Noise` instance is passes as first argument, its autocovarince is
estimated and passed to the fitting function.

Note: Change initial guess if fit does not converge to sensible result.
"""
function fit_noise(
    lags::Vector{T},
    autocov::Vector{T};
    α_guess=1.0,
    λ_guess=1.0
) where {T<:AbstractFloat}
    fit = curve_fit(gauss_kernel, lags, autocov, [α_guess, λ_guess])
    return MvGaussianNoiseModel{T}(lags[1], fit.param...)
end

function fit_noise(n::Noise{T}; α_guess=1.0, λ_guess=1.0) where {T}
    lags, acov = estimate_autocov(n)
    fit_noise(lags, acov; α_guess=α_guess, λ_guess=λ_guess)
end