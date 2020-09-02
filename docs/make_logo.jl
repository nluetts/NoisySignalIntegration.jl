using Colors
using Plots
using NoisySignalIntegration
using Random: seed!

gausspeaks(x, p) = sum([@. A * 1/√(2π * σ^2) * exp(-(x - μ)^2 / (2σ^2)) for (A, μ, σ) in p])

c1, c2, c3 = let
    seed!(1)
    x = 0:0.1:100
    y = gausspeaks(x, [(3.3, 30.0, 4.0), (4.1, 70.0, 4.0), (1, 50.0, 20.0)])
    c = Curve(x, y)
    uc = add_noise(c, MvGaussianNoiseModel(0.1, 0.01, 1.0), 3)
    c1, c2, c3 = [NoisySignalIntegration.get_draw(i, uc) for i in 1:3]
    c2 += 0.5
    c3 += 1
    c1, c2, c3
end

p = let
    r = [203, 60, 51]  ./ 255
    v = [149, 88, 178] ./ 255
    b = [64, 99, 216]  ./ 255
    g = [56, 152, 38]  ./ 255
    plot(; background_color = :transparent, color=:black, ticks=:none, legend=false, linewidth=3, axis=false);
    plot!(c1, 18., 40., fillcolor=RGB(r...), fillalpha=1);
    plot!(c1, 60., 84., fillcolor=RGB(v...), fillalpha=1);
    plot!(c2, 19., 42., fillcolor=RGB(b...), fillalpha=1);
    plot!(c2, 59., 84., fillcolor=RGB(g...), fillalpha=1);
    plot!(c1; color=:black, linewidth=8);
    plot!(c2; color=:black, linewidth=8)
    plot!(annotations=(70, 0.13, Plots.text("?", 50, RGBA(1, 1, 1, 0.25))))
end

savefig(p, joinpath(@__DIR__, "src", "assets", "logo.png"))