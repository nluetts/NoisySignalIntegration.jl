using LsqFit: curve_fit
using Polynomials: fit
using StatsBase: autocov

"""
    detrend(x, y, poly_order)

Subtract polynomial baseline from y data.
"""
detrend(x, y, poly_order) = y - fit(x, y, poly_order)(x)

"""
    function estimate_noise_parameters(slc::Slice; p₀=[1.0, 1.0], tolG=1e-20,
                                       detrend_order=0)

Estimate parameters α and λ of the exponentiated quadratic covariance function
k(Δx) = α² exp(-0.5 (Δx/λ)²) by fitting the estimated autocovaraiance of the
noise n. Remove a slow varying background of the noise by subtracting a
polynomial of order `detrend_order` if detrend_order > 0.
"""
function estimate_noise_parameters(slc::Slice; p₀=[1.0, 1.0], g_tol=1e-20,
                                  detrend_order=0)
    x, n = begin
        spec = crop(slc)
        spec.x, spec.y
    end
    @. k(Δx, β) = (β[1]^2)*exp(-0.5*(Δx/β[2])^2)
    δx = abs(x[1]-x[2]) # x-values must be evenly spaced
    lags = 1:50
    n_ = detrend_order > 0 ? detrend(x, n, detrend_order) : n
    acov = autocov(n_, lags)
    Δx = collect(δx*lags)

    fit = curve_fit(k, Δx, acov, p₀, g_tol=g_tol)

    fig = plot(xlabel="lag", ylabel="autocovariance")
    plot!(fig, Δx, acov, label="estimate")
    plot!(fig, Δx, k(Δx, fit.param), label="fit")

    fig, fit.param
end