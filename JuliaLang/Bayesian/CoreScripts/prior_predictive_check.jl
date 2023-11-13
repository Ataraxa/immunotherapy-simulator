#= File Description

This file contains a script to run a prior predictive check on a model.
User must specify both the priors and likelihood.
=#

using Plots: plot, plot!
using Turing
using StatsBase: percentile, mode, median
using Distributions


include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_restricted.jl")

### Script Settings
num_samples = 1_000
simul_matrix = Matrix{Float64}(undef, num_samples, 271) # (undef, row, col)
problem = create_problem()
### Priors 
ln_k₆_prior = truncated(Cauchy(0, 1); lower=-100, upper=0)
ln_d₁_prior = truncated(Cauchy(0, 1); lower=0, upper=7)
ln_s₂_prior = truncated(Cauchy(0, 1); lower=-100, upper=0)


### Main
ln_k₆ = rand(ln_k₆_prior, num_samples)
ln_d₁ = rand(ln_d₁_prior, num_samples)
ln_s₂ = rand(ln_s₂_prior, num_samples)

for i in 1:num_samples
    p = [ln_k₆[i], ln_d₁[i], ln_s₂[i]] .|> exp 

    re_p, u0 = repack_params(updateParams3(p...))
    predictions = solve(problem; p=re_p, saveat=0.1)
    simul_matrix[i,:] = predictions[4,:] + predictions[5,:]
end

plot(0:0.1:27.0, simul_matrix'; legend=false, color="#BBBBBB", alpha=0.3)
plot!(
    0:0.1:27, percentile.(eachcol(simul_matrix), 5);
    fillrange =  percentile.(eachcol(simul_matrix), 95),
    fillalpha = 0.3,
    alpha = 0
    )
plot!(0:0.1:27.0, median.(eachcol(simul_matrix)))

    