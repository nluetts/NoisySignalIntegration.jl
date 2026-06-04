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
    crv + randn(length(crv)) .* maximum(crv.y) .* 0.05
end

function main()
    # Load the test dataset
    dataset = NoisySignalIntegration.testdata_1()
    udataset = add_noise(dataset, MvGaussianNoiseModel(0.1, 0.05, 0.5)) # make sure this fits testdata_1()

    bds = [
        # TODO: for some reason, fitting the first bound/peak takes
        # a lot more iterations and time than the second. Why is that?
        UncertainBound(Uniform(11.0, 12.0), Uniform(18.0, 19.0)),
        UncertainBound(Uniform(24.0, 26.0), Uniform(34.0, 36.0)),
    ]

    res = mc_fit(udataset, bds)

    # Plot the dataset along with both fitted curves
    udataset, res, bds
end
