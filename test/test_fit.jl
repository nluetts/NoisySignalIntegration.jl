using NoisySignalIntegration
using NoisySignalIntegration: Curve, UncertainCurve, UncertainBound, fwhm
using Distributions: Uniform
using Test
using Random
using Plots

# Widths of test spectrum
W1 = 2.0
W2 = 1.0
W3 = 1.5

# Gaussian functions defined by height and fwhm
function gauss(x, x0, A, fwhm)
    sigma = fwhm / sqrt(8*log(2))
    @. A * exp(-(x - x0)^2 / (2 * sigma^2))
end

# Draw three Gaussians
function test_spectrum()
    xs = collect(0:0.025:100)
    ys = zeros(length(xs))
    ys += gauss(xs, 20, 1, W1)
    ys += gauss(xs, 50, 2.5, W2)
    ys += gauss(xs, 70, 1.5, W3)
    Curve(xs, ys)
end

# ... plus noise
function test_spectrum_noisy()
    crv = test_spectrum()
    Random.seed!(42)
    crv + randn(length(crv)) .* maximum(crv.y) .* 0.01
end

function main()
    @testset "regression tests" begin
        @testset "fit of single peak" begin
            c = test_spectrum_noisy()
            uc = add_noise(c, GaussianNoiseModel(0.01))
            ub = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 15.0, 25.0), uc)
            res = mc_fit(uc, ub)
            return res[1], c, uc, ub
        end
    end
end
