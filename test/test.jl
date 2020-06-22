##
using Revise

##
using MCIntegrate


## Imports
using Distributions: MvNormal
using Plots
using Random

## Test spectrum

curve = Curve([1.0, 2.0], [3.0, 4.0])
@assert length(curve) == 2
@assert curve[end] == (2.0, 4.0)
@assert [curve; curve] == Curve([1.0, 2.0, 1.0, 2.0], [3.0, 4.0, 3.0, 4.0])
@assert curve + curve.y == Curve([1.0, 2.0], [6.0, 8.0])

curve = Noise([1.0, 2.0], [4.0, 4.0])
@assert length(curve) == 2
@assert curve[end][1] == 2.0
@assert curve[end][2] + 1 ≈ 1
@assert [curve; curve] == Noise([1.0, 2.0, 1.0, 2.0], [4.0, 4.0, 4.0, 4.0])

##

function get_test_spectrum(seed)
    Random.seed!(seed)
    x = collect(0:0.1:100)
    y = @. exp(-(x-15)^2) + 2 * exp(-(x-30)^2)
    @. y += 1.0 + x*1.5e-2 - (x-50)^3*3e-6
    n = length(x)
    δx = x[2] - x[1] 
    return Curve(x, y .+ (get_cov(δx, n, 0.1, 0.5) |> MvNormal |> rand))
end

spec = get_test_spectrum(1)
plot(spec)

## cut out parts for analysis

slc_bands = crop(spec, 5.0, 40.0)
slc_noise = crop(spec, 40.0, 100.0)

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
b1 = WidthBound(15.0, scale_shift_beta(2.0, 2.0, 6.0, 7.0))
b2 = LeftRightBound(scale_shift_beta(2.0, 2.0, 25.0, 28.0), scale_shift_beta(2.0, 2.0, 35.0, 38.0))

p = histogram(sample(b2, 100000))

##
integral_samples = mc_integrate(spec, noise_param, [b1, b2]; N=10_000)  # left off here
