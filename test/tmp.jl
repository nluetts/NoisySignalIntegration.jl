module Tmp

using Plots
using MCIntegrate
using Distributions: Beta, Normal

d1 = Beta(2, 2)
d2 = Normal(0, 1)
md1 = ScaledShiftedCUD(d1, 3.0, 0.1, 0.5)
md2 = ScaledShiftedCUD(d2, 4.0, 0.1, 0.0)
ub = UncertainBounds(md1, md2)

p = histogram(sample(ub, 100000))

end