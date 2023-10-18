using Plots: plot, plot!
using Distributions
using StatsPlots
using ForwardDiff: Dual

include("../Model/Differential/dde_to_bayesian.jl")
include("../Library/optimisation_lib.jl")

pb = restricted_dde_space()
init_k6 = 0.4524371106042169
init_d1 = 17.493801636164296
init_s2 = 0.27859689121549636

fd_k6 = Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}}(0.4524371106042169, 1.651151666195318, 0.0, 0.0)
fd_d1 = Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}}(17.493801636164296, 0.0, 2.339124341340631, 0.0)
fd_s2 = Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}}(0.2785968912154963, 0.0, 0.0, 1.141240725012459)
p = [fd_k6, fd_d1, fd_s2]
sol = solve(pb; p=p, saveat=0.1)
pred_vol = sol[4,:] + sol[5,:]
selected_days = [7,10,12,14,16,18,20,25,27]
s = 0.1
sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]
plot(sol.t, [sol[4,i].value + sol[5,i].value for i in eachindex(sol)])

p = [fd_k6.value, fd_d1.value, fd_s2.value]
sol = solve(pb; p=p, saveat=0.1)
pred_vol = sol[4,:] + sol[5,:]
selected_days = [7,10,12,14,16,18,20,25,27]
s = 0.1
sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]
plot!(sol.t, sol[4,:] + sol[5,:])


# plot!(selected_days, sliced_pred)
# plot!(cbd_il_12["active_days"], cbd_il_12["mean"])


# normal = truncated(Normal(11, 3); lower=0, upper=20)
# plot(density(rand(normal, 1000)))

# logi = truncated(Normal(log(11), log(1)); lower=log(0), upper=log(20))
# plot!(density(exp.(rand(logi, 1000))))
