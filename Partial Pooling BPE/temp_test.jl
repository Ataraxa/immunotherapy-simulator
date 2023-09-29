using Plots: plot, plot!

include("./Model/Differential/dde_to_bayesian.jl")
include("./Library/optimisation_lib.jl")

pb = restricted_dde_space()
k6 = 0.03113
d1 = 3.25402
s2 = 7.18949
p = [exp(k6), exp(d1), exp(s2)]
print(p)
sol = solve(pb; p=p, saveat=0.1)
plot(sol.t, sol[4,:]+sol[5,:])
# plot!(cbd_il_12["active_days"], cbd_il_12["mean"])
