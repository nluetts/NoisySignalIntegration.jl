@recipe function plot_recipe(crv::AbstractCurve)
    xguide --> "x"
    yguide --> "y"
    return crv.x, crv.y
end