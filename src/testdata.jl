gausspeaks(x, p) = sum([@. A * 1/√(2π * σ^2) * exp(-(x - μ)^2 / (2σ^2)) for (A, μ, σ) in p])

"""
    generate_testdata(grid, peaks, nm::AbstractNoiseModel; seedvalue=1) :: Curve

Return simulated test data (sum of Gauss peaks plus noise).

# Arguments

- `grid`: grid on which simulation shall be calculated
- `peaks`: list of tuples, each holding area, position, and width (standard deviation) of Gauss peaks
- `nm :: AbstractNoiseModel`: model that controls noise in simulation
"""
function generate_testdata(grid, peaks, nm::AbstractNoiseModel; seedvalue=1, baseline=nothing)
    seed!(seedvalue)
    # simulate FTIR spectrum with two bands
    # the true intensity ratio is 1 : 2
    signals = gausspeaks(grid, peaks)
    baseline = baseline === nothing ? zeros(eltype(grid), length(grid)) : baseline
    uc = add_noise(Curve(grid, baseline + signals), nm, 1)
    return get_draw(1, uc)
end


"""
Return test dataset with two peaks with area ratio 1 : 2 and correlated noise
(α = 0.05, λ = 0.5).
"""
testdata_1() = let
    x = collect(0:0.1:100)
    baseline = @. 1.0 + x*1.5e-2 - (x-50)^3*3e-6
    signals = [(1, 15, √0.5), (2, 30, √0.5)]
    generate_testdata(x, signals, MvGaussianNoiseModel(0.1, 0.05, 0.5); baseline=baseline)
end


"""
Return test dataset with two peaks with area ratio 1 : 2 and uncorrelated noise (σ = 0.1).
"""
testdata_2() = let
    x = collect(0:0.1:100)
    baseline = @. 1.0 + 1.5e-4x^2 - 3e-6x^3
    signals = [(1, 15, √0.5), (2, 30, √0.5)]
    generate_testdata(x, signals, GaussianNoiseModel(0.1); baseline=baseline)
end


"""
Return test dataset from testdata_1() with uneven grid spacing.
"""
testdata_3() = let
    x = collect(0:0.1:100)
    x = @. x + 1e-3x^2
    baseline = @. 1.0 + x*1.5e-2 - (x-50)^3*3e-6
    signals = [(1, 15, √0.5), (2, 30, √0.5)]
    generate_testdata(x, signals, MvGaussianNoiseModel(0.1, 0.05, 0.5); baseline=baseline)
end


"""
Return test dataset with four peaks with area ratio 1 : 2 : 3 : 0.5 and correlated noise
(α = 0.05, λ = 0.5).
"""
testdata_4() = let
    x = collect(0:0.1:200)
    baseline = @. 1.0 + x*2.5e-3 - (x-50)^3*8e-8
    signals = [(1, 15, √0.5), (2, 30, √0.5), (3, 60, √0.5), (1, 85, √0.5)]
    generate_testdata(x, signals, MvGaussianNoiseModel(0.1, 0.05, 0.5); baseline=baseline, seedvalue=42)
end



