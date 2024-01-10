# Defines the UncertainBound type to handle uncertain untegration start and end points


# Define which quantiles are stored in UncertainBound objects
const BOUND_QUANTILES = 0.1:0.2:0.9 |> collect
const IDX_BOUND_MEDIAN = 3 # the median is at index 3


"""
    UncertainBound{T, N}

Used to define uncertain integration bounds.

Holds `N` samples of start and end points of type `T` for numeric integration.

# Fields

- `left::Particles{T, N}`: samples of starting point
- `right::Particles{T, N}`: samples of end point

# Constructors

    UncertainBound(left::S, right::T, N::Int=10_000) where {S <: ContinuousUnivariateDistribution, T <: ContinuousUnivariateDistribution}

Create an `UncertainBound` from two distributions `left` and `right` defining uncertain start and end points of the integration window.

**Example**

Start point falls in the range [1, 2] with uniform probability, end point falls in the range [5, 6] with uniform probability:

```jldoctest UBexample
julia> using Distributions

julia> using Random: seed!; seed!(1);

julia> UncertainBound(Uniform(1, 2), Uniform(5, 6), 20_000)
UncertainBound{Float64, 20000}(start = 1.5 ± 0.29, end = 5.5 ± 0.29)
```

---

    UncertainBound(pos::T, width::ContinuousUnivariateDistribution, uc::UncertainCurve{T, N}) where {T, N}

Create an `UncertainBound` from a position `pos`, a distribution for the integration window width `width`, and an uncertain curve object `uc`.
Useful for integrating symmetric peaks.
The constructor will find start and end points for the integration window such that, for each sample in `uc`, the integration window will
be symmetric around the peak that falls in the interval `pos` ± (sample of `width`)/2, i.e. the integration window will "follow" the peak
position from draw to draw (note that, due to noise, the peak may change its position from draw to draw).

**Example**

Uncertainty of width is modeled with a scaled and shifted Beta(2, 2) distribution:

```jldoctest UBexample
julia> uc = begin # create uncertain curve with one symmetric peak
           x = 0:0.1:10;
           y = @. exp(-(x - 5)^2)
           add_noise(Curve(x, y), MvGaussianNoiseModel(0.1, 0.03, 0.5))
       end;

julia> ub = UncertainBound(5., scale_shift_beta(2, 2, 3.5, 4.0), uc)
UncertainBound{Float64, 10000}(start = 3.12 ± 0.063, end = 6.87 ± 0.064)
```

---

    UncertainBound(pos::Vector{T}, width::ContinuousUnivariateDistribution, uc::UncertainCurve{T, N}) where {T, N}

Create several `UncertainBound` objects with identical integration window widths in each draw.
Useful for integration of several symmetric peaks for which the same width may be assumed.


**Example**

Uncertainty of width is modeled with a scaled and shifted Beta(2, 2) distribution:

```jldoctest UBexample
julia> uc = begin # create uncertain curve with one symmetric peak
           x = 0:0.1:10;
           y = @. exp(-(x - 3)^2/0.15) + exp(-(x - 7)^2/0.15)
           add_noise(Curve(x, y), MvGaussianNoiseModel(0.1, 0.03, 0.5))
       end;

julia> ubs = UncertainBound([3., 7.], scale_shift_beta(2, 2, 1.3, 1.5), uc)
2-element Vector{UncertainBound{Float64, 10000}}:
 UncertainBound{Float64, 10000}(start = 2.3 ± 0.022, end = 3.7 ± 0.022)
 UncertainBound{Float64, 10000}(start = 6.3 ± 0.022, end = 7.7 ± 0.022)
```
"""
struct UncertainBound{T,N}
    left::Particles{T,N}
    right::Particles{T,N}
    # quantiles for local baseline subtraction and plotting
    _left_quantiles
    _right_quantiles
end

function Base.show(io::IO, c::UncertainBound{T,N}) where {T,N}
    print(io, "UncertainBound{$T, $N}(start = $(c.left), end = $(c.right))")
end


# Constructors

# Base constructor that injects quantiles
function UncertainBound(left::Particles{T,N}, right::Particles{T,N}) where {T,N}
    qs = BOUND_QUANTILES
    lqs, rqs = [[pquantile(lr, q) for q in qs] for lr in (left, right)]
    return UncertainBound{T,N}(left, right, lqs, rqs)
end


# Create a left/right bound
function UncertainBound(left::S, right::T, N::Int=10_000) where {S<:ContinuousUnivariateDistribution,T<:ContinuousUnivariateDistribution}
    left = Particles(N, left)
    right = Particles(N, right)
    return UncertainBound(left, right)
end


# Create several, correlated width bounds
function UncertainBound(
    pos::Vector{T},
    width::ContinuousUnivariateDistribution,
    uc::UncertainCurve{T,N}
) where {T,N}

    M = length(pos)
    left = Array{T}(undef, M, N)
    right = Array{T}(undef, M, N)
    width = Particles(N, width)

    for i ∈ 1:M
        pᵢ = pos[i]
        for j ∈ 1:N
            cⱼ = get_draw(j, uc) # sample j from curve
            wⱼ = get_draw(j, width) # sample j from width
            wⱼ < zero(T) && replace("""Negative integration window width encountered ($wⱼ),
                                       make sure the support region of your width distribution
                                       does not allow negative values!""", "\n" => " ") |> ArgumentError |> throw
            left[i, j], right[i, j] = left_right_from_peak(uc.x, cⱼ.y, pᵢ, wⱼ)
        end
    end

    return [UncertainBound(Particles(left[i, :]), Particles(right[i, :])) for i in 1:M]
end


# Create single width bound
function UncertainBound(
    pos::T,
    width::ContinuousUnivariateDistribution,
    uc::UncertainCurve{T,N}
) where {T,N}
    bnd = UncertainBound([pos], width, uc)
    return bnd[1]
end

"""
    get_draw(n, bnd::UncertainBound)

Retrieve the `n`th sample of the samples stored in `UncertainBound` `bnd`.
"""
get_draw(n, bnd::UncertainBound) = [get_draw(n, bnd.left), get_draw(n, bnd.right)]

"""
    mean(uc::UncertainBound)

Retrieve the mean of the `UncertainBound` start and end point.
"""
Statistics.mean(bnd::UncertainBound) = [Statistics.mean(bnd.left), Statistics.mean(bnd.right)]

"""
    scale_shift_beta(α, β, a, b)

Create a scaled and shifted `Beta(α, β)` distribution.
Samples fall in the interval [`a`, `b`].
"""
function scale_shift_beta(α, β, a, b)
    return Beta(α, β) * (b - a) + a
end
