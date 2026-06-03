using NoisySignalIntegration
using NoisySignalIntegration: Curve, UncertainCurve, UncertainBound, fwhm
using Test
using Random

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


@testset "fwhm of simple gaussian curves" begin
    @info "Testing FWHM determination of simple Gaussians ..."
    crv = test_spectrum()
    @test fwhm(crv, 10., 30.).full_width == W1
    @test fwhm(crv, 40., 60.).full_width == W2
    @test fwhm(crv, 60., 80.).full_width == W3
end

@testset "fwhm of simple gaussian curves" begin
    @info "Testing FWHM determination of simple Gaussians on slope ..."
    crv = test_spectrum()
    crv += crv.x * 0.001
    @test fwhm(crv, 10., 30.; local_baseline=true).full_width == W1
    @test fwhm(crv, 40., 60.; local_baseline=true).full_width == W2
    @test fwhm(crv, 60., 80.; local_baseline=true).full_width == W3
end
