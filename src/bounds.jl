"""
   draw_from_beta(draws, μ, s)

Return vector of random draws from a beta(2,2) distribution.

Return values are shifted by μ and scaled by s and fall in the interval
[μ - ½s, μ + ½s].
"""
function draw_from_beta(draws, μ, s) :: Vector{Float64}
    (rand(Beta(2.0, 2.0), draws) .- 0.5) .* s .+ μ
end