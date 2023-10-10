using Plots: plot, plot!
using Distributions
using StatsPlots
using ForwardDiff: Dual

include("./Model/Differential/dde_to_bayesian.jl")
include("./Library/optimisation_lib.jl")

pb = restricted_dde_space()
k6 = Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}}(3.0633312778684614,5.642811918390731,0.0,0.0).value
d1 = Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}}(10.097704850342026,0.0,6.855223618900589,0.0).value
s2 = Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}}(0.01849635637941573,0.0,0.0,0.1204497285489698).value
p = [k6, d1, s2]
println(p)
sol = solve(pb; p=p, saveat=0.1)
pred_vol = sol[4,:] + sol[5,:]
selected_days = [7,10,12,14,16,18,20,25,27]
s = 0.1
sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]
println(sliced_pred)
plot(sol.t, sol[4,:] + sol[5,:])
# plot(sol.t, [sol[4,i].value + sol[5,i].value for i in eachindex(sol)])
# plot!(selected_days, sliced_pred)
# plot!(cbd_il_12["active_days"], cbd_il_12["mean"])


# normal = truncated(Normal(11, 3); lower=0, upper=20)
# plot(density(rand(normal, 1000)))

# logi = truncated(Normal(log(11), log(1)); lower=log(0), upper=log(20))
# plot!(density(exp.(rand(logi, 1000))))
