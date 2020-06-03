##
using Revise

##
using MCIntegrate


## Imports
using Distributions
using Plots
using Random

## Test spectrum

spec = Spectrum([1.0, 2.0], [3.0, 4.0])
@assert length(spec) == 2
@assert spec[end] == (2.0, 4.0)
@assert [spec; spec] == Spectrum([1.0, 2.0, 1.0, 2.0], [3.0, 4.0, 3.0, 4.0])

##

function get_test_spectrum(seed)
    Random.seed!(seed)
    x = collect(0:0.1:100)
    y = @. exp(-(x-15)^2) + 2 * exp(-(x-30)^2)
    @. y += 1.0 + x*1.5e-2 - (x-50)^3*3e-6
    return Spectrum(x, y .+ (get_cov(x, 0.1, 0.5) |> MvNormal |> rand))
end

spec = get_test_spectrum(1)
plot(spec)

## cut out parts for analysis

slc_bands = crop(spec, 5.0, 40.0)
slc_noise = crop(spec, 40.0, 100.0)

plot(slc_bands, alpha=0.5) |> (x -> plot!(x, slc_noise, alpha=0.5))

##

noise_sample = NoiseSample(slc_noise)