using Plots: plot, plot!
using Distributions
using StatsPlots

include("./Model/Differential/dde_to_bayesian.jl")
include("./Library/optimisation_lib.jl")

# pb = restricted_dde_space()
# k6 = 0.03113
# d1 = 3.25402
# s2 = 7.18949
# p = [exp(k6), exp(d1), exp(s2)]
# print(p)
# sol = solve(pb; p=p, saveat=0.1)
# plot(sol.t, sol[4,:]+sol[5,:])
# # plot!(cbd_il_12["active_days"], cbd_il_12["mean"])


normal = truncated(Normal(11, 3); lower=0, upper=20)
plot(density(rand(normal, 1000)))

logi = truncated(Normal(log(11), log(1)); lower=log(0), upper=log(20))
plot!(density(exp.(rand(logi, 1000))))
