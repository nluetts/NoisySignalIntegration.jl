using Revise

##
using MCIntegrate


## Imports
using Distributions: MvNormal
using Plots
using Random
using MonteCarloMeasurements

const SAMPLES = 50_000

function get_test_spectrum(seed)
    Random.seed!(seed)
    x = collect(0:0.1:200)
    y = @. exp(-(x-15)^2) + 2 * exp(-(x-30)^2)
    @. y += 1.0 + x*1.5e-2 - (x-50)^3*3e-7
    n = length(x)
    δx = x[2] - x[1] 
    return Curve(x, y .+ (get_cov(δx, n, 0.1, 0.5) |> MvNormal |> rand))
end

spec = get_test_spectrum(1);
noise = crop(spec, 50, 200)
spec = crop(spec, 10, 40)

unc_spec = UncertainCurve(spec.x, spec.y + correlated_noise(noise, length(spec), SAMPLES));

ubnd1 = UncertainBound(
    Particles(SAMPLES, scale_shift_beta(2.0, 2.0, 13.0, 13.5)),
    Particles(SAMPLES, scale_shift_beta(2.0, 2.0, 16.5, 17.0))
)
ubnd2 = UncertainBound(
    Particles(SAMPLES, scale_shift_beta(2.0, 2.0, 28.0, 28.5)),
    Particles(SAMPLES, scale_shift_beta(2.0, 2.0, 31.5, 32.0))
)

ubnd3 = UncertainBound(
    15.0,
    scale_shift_beta(2.0, 2.0, 3.0, 4.0),
    unc_spec
)

ubnd4, ubnd5 = UncertainBound(
    [15.0, 30.0],
    scale_shift_beta(2.0, 2.0, 3.0, 4.0),
    unc_spec
)

_get_N_samples(::Particles{T, N}) where {T, N} = N

areas = mc_integrate(unc_spec, [ubnd1, ubnd2, ubnd3, ubnd4, ubnd5])
