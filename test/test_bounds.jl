module Tmp

using Plots
using MCIntegrate
using Distributions: Beta, Normal

md1 = scale_shift_beta(2.0, 2.0, 3.0, 6.0)
md2 = Normal(0.0, 0.1)
ub = UncertainBounds(md1, md2)

p = histogram(sample(ub, 100000))

end