using .Plots

"""
    animate_draws(uc::UncertainCurve, bnds::Vector{UncertainBound}; n=20, fps=5, filepath=nothing, kw...)

Create a gif animation showing the first `n` random draws.
Requires the Plots.jl package.

Create an animation of the Monte-Carlo iterations when integrating the
`UncertainCurve` `uc` using the `UncertainBound`s `bnds`.

# Keyword arguments

`n`: number of draws that are included in animation

`fps`: frames-per-second

`filepath`: if provided, the gif-animation will be stored using this filepath,
otherwise a temporary file will be created.

Further keyword arguments are passed on the `Plots.plot()` function.
"""
function animate_draws(
    uc::T, bnds::Vector{S}; n=20, fps=5, filepath=nothing, label="", kw...
) where {T<:UncertainCurve, S<:UncertainBound}
    # get minimum and maximum y value of all curve draws
    # in the selected `n` samples, used to determine the
    # y-span when plotting
    miny, maxy = [
        mm(curve.y for curve ∈ [NoisySignalIntegration.get_draw(i, uc) for i ∈ 1:n]) |> mm
        for mm in (minimum, maximum)
    ]
    anim = @animate for i in 1:n
        plot(NoisySignalIntegration.get_draw(i, uc), bnds, i; label=label, ylim=(miny*0.9, maxy*1.1), kw...)
    end
    fp = filepath === nothing ? tempname()*".gif" : filepath
    return gif(anim, fp, fps=fps)
end

export animate_draws