# --------------------------------------------
# enable plotting of custom types via recipes
# --------------------------------------------


# --------------------------------------------
# enable plotting of curves and noise samples
# --------------------------------------------
@recipe function plot_recipe(crv::AbstractCurve)
    xguide --> "x"
    yguide --> "y"
    return crv.x, crv.y
end