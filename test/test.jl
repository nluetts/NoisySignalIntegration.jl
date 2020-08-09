"""
This test file can be run manually to check some of the
plotting functionallity.
"""

##
using Revise

##
using MCIntegrate


## Imports
using Distributions: MvNormal
using Plots
using Random


function get_test_spectrum(seed)
    Random.seed!(seed)
    x = collect(0:0.1:200)
    y = @. exp(-(x-15)^2) + 2 * exp(-(x-180)^2)
    @. y += 1.0 + x*1.5e-2 - (x-50)^3*3e-7
    n = length(x)
    δx = x[2] - x[1] 
    return Curve(x, y .+ (get_cov(δx, n, 0.1, 0.5) |> MvNormal |> rand))
end

spec = get_test_spectrum(1)
plot(spec)

## cut out parts for analysis

slc_band1 = crop(spec, 10.0, 20.0)
slc_band2 = crop(spec, 175.0, 185.0)
slc_bands = [slc_band1; slc_band2]
slc_noise = crop(spec, 25.0, 170.0)

plot(slc_bands, alpha=0.5) |> (x -> plot!(x, slc_noise, alpha=0.5))

## fit noise and plot

noise_sample = Noise(slc_noise, 3)
noise_param = fit_noise(noise_sample)

plot_autocov(noise_sample, noise_param)

## 
plot(noise_param; noise_samples=3)

##
Random.seed!(42)
plot(noise_sample, noise_param)

##
b1 = WidthBound(180.0, scale_shift_beta(2.0, 2.0, 3.0, 4.0), SUBTRACT_LOCAL)
b2 = LeftRightBound(scale_shift_beta(2.0, 2.0, 13.0, 13.5), scale_shift_beta(2.0, 2.0, 16.5, 17.0), INCLUDE)
b3 = LeftRightBound(scale_shift_beta(2.0, 2.0, 178.0, 178.5), scale_shift_beta(2.0, 2.0, 181.5, 182.0), INCLUDE)
b4 = clone(b1, 15.0)

p = histogram(sample(b2, 100000))

##

plot(slc_bands, [b1, b2, b3, b4], noise_param)

##
integral_samples = mc_integrate(slc_bands, [b1, b2, b3, b4], noise_param; N=100_000)

## 
histogram(integral_samples, alpha=0.5, normalize=true)

histogram(integral_samples[:, 1]./integral_samples[:, 4])
histogram!(integral_samples[:, 3]./integral_samples[:, 2], alpha=0.5)

## 
x = collect(1:1.0:10)
y = [-1.5, 0.0, 1.0, 2.0, 2.5, 2.0, 2.0, 1.0, 0.0, -1.5]
c = Curve(x, y)
scatter(c)
plot!(c, 0.0, 9.5)