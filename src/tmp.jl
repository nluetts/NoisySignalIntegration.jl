#######################################################
### plot integration bounds ###########################
#######################################################

function _plot_bound!(p::Plots.Plot, crv::Curve{T}, left::T, right::T) where {T}
    bounds_parameters = get_integration_bounds(crv.x, crv.y, left, right)
    l, r = bounds_parameters[2:3]
    x_l, x_r, y_l, y_r = bounds_parameters[12:15]
    x_ = [x_l; crv.x[l+1:r-1]; x_r; x_l]
    y_ = [y_l; crv.y[l+1:r-1]; y_r; y_l]
    return plot!(p, x_, y_; fill=(0, 0.5, :orange), linewidth=0, label=nothing)
end

function Plots.plot!(
    p::Plots.Plot,
    crvs::Vector{T},
    bnds::Vector{S}
) where {T <: AbstractCurve, S <: AbstractUncertainBound}
    for c in crvs
        plot!(p, c, label=nothing)
        for b in bnds
            l, r = typeof(b) <: WidthBound ? sample(b, c) : sample(b)
            _plot_bound!(p, c, l, r)
        end
    end
    return p
end
Plots.plot(crvs::Vector{AbstractCurve}, bnds::Vector{AbstractUncertainBound}) = Plots.plot!(plot(), crvs, bnds)
Plots.plot!(p::Plots.Plot, crv::AbstractCurve, bnds::Vector{AbstractUncertainBound}) = Plots.plot!(p, [crv], bnds)
Plots.plot(crv::AbstractCurve, bnds::Vector{AbstractUncertainBound}) = Plots.plot!(plot(), [crv], bnds)
Plots.plot!(p::Plots.Plot, crvs::Vector{AbstractCurve}, bnd::AbstractUncertainBound) = Plots.plot(p, crvs, [bnd])
Plots.plot(crvs::Vector{AbstractCurve}, bnd::AbstractUncertainBound) = Plots.plot(plot(), crvs, [bnd])
Plots.plot!(p::Plots.Plot, crv::AbstractCurve, bnd::AbstractUncertainBound) = Plots.plot!(p, [crv], [bnd])
Plots.plot(crv::T, bnd::S) where {T <: AbstractCurve, S <: AbstractUncertainBound} = Plots.plot(plot(), [crv], [bnd])