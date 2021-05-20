# Baseline Handling

## Build-in

- none
- local_baseline
- subtract_baseline

```@eval
using NoisySignalIntegration
using Plots
using Random: seed!

seed!(1)
nsi = NoisySignalIntegration

let
    n = 20
    c = crop(nsi.testdata_1(), 0, 50)
    uc = add_noise(c, GaussianNoiseModel(0.1))
    ubleft = UncertainBound(15.0, scale_shift_beta(2.0, 2.0, 4.0, 5.0), uc)
    ubright = UncertainBound(30.0, scale_shift_beta(2.0, 2.0, 6.0, 7.0), uc)

    spany = [
        mm(curve.y for curve ∈ [nsi.get_draw(i, uc) for i ∈ 1:n]) |> mm
        for mm in (minimum, maximum)
    ]
    spany = (spany[1]*0.9, spany[2]*1.1)
    anim = @animate for i in 1:n
        kw = Dict(:ylim => spany, :legend => :topleft)
        p1 = plot(nsi.get_draw(i, uc), [ubleft, ubright], i; subtract_baseline=true, label=["subtract_baseline" "" ""], kw...)
        p2 = plot(nsi.get_draw(i, uc), [ubleft, ubright], i; local_baseline=true, label=["local_baseline" "" ""], kw...)
        plot(p1, p2; layout=(2, 1), size=(400, 300))
    end
    gif(anim, "baseline_anim.gif", fps=5)
    nothing
end
```

![build-in baseline handling animation](baseline_anim.gif)

## Preprocessing

- without uncertainty
- with uncertainty