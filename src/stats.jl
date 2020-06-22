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

# GaussianNoiseModel

function sample!(s::AbstractArray{T}, nm::GaussianNoiseModel{T}) where {T}
    randn!(s)
    s *= nm.σ
    return s
end

function sample(nm::GaussianNoiseModel{T}, m::Integer=100, n::Integer=1) where {T}
    s = Array{T}(undef, m, n)
    return sample!(s, nm)
end

# MvGaussianNoiseModel

function sample!(
    s::AbstractArray{T},
    nm::MvGaussianNoiseModel{T},
) where {T<:AbstractFloat}
    m = size(s)[1]
    rand!(MvNormal(get_cov(nm.δx, m, nm.α, nm.λ)), s)
    return s
end
sample(nm::MvGaussianNoiseModel{T}, m::Integer=100, n::Integer=1) where {T <: AbstractFloat}= sample!(Array{T}(undef, m, n), nm)

# LeftRightBound

function sample!(s::AbstractArray{T}, b::LeftRightBound) where {T <: AbstractFloat}
    for i in eachindex(s[:, 1])
        s[i, 1] = rand(b.left)
        s[i, 2] = rand(b.right)
    end
    return s
end

function sample(b::LeftRightBound, samples::Integer=1)
    s = Array{Float64}(undef, samples, 2)
    return sample!(s, b)
end

# WidthBound

function sample!(spls::AbstractArray{T}, wb::WidthBound{T}, s::Curve{T}) where {T <: AbstractFloat}
    rand!(wb.width, spls[:, 1])
    n, _ = size(spls)
    for i in 1:n
        spls[i, :] = left_right_from_peak(s.x, s.y, wb.loc, spls[i, 1])
    end
    return spls
end

function sample(wb::WidthBound{T}, s::Curve{T}, samples::Integer=1) where {T <: AbstractFloat}
    spls = Array{T}(undef, samples, 2)
    return sample!(spls, wb, s)
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