## Test to implement a try loop to avoid maxiters errors

include("../Model/mechanistic_model.jl")
include("../Model/Bayesian/priors.jl")

pb = create_problem()
p = christian_true_params
var_p = [11, 12, 21]
p[var_p] .= exp.([7, -5, 7])
println(p)

pred = solve(pb; p=p, saveat=0.1)