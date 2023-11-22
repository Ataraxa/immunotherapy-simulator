#= File Description

This file contains a script to run a prior predictive check on a model.
User must specify both the priors and likelihood.
=#

using Plots: plot, plot!
using Turing
using StatsBase: percentile, mode, median
using Distributions
using DelimitedFiles

include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_restricted.jl")

### Script Settings
num_samples = 500
simul_matrix = Matrix{Float64}(undef, num_samples, 271) # (undef, row, col)
problem = create_problem()

### Priors 
ln_k₆_prior = truncated(Cauchy(-0.61, 1); lower=-100, upper=0)
ln_d₁_prior = truncated(Cauchy(3, 1); lower=0, upper=7)
ln_s₂_prior = truncated(Cauchy(-1.0, 1); lower=-100, upper=0)

### Main
ln₍k₆₎ = rand(ln_k₆_prior, num_samples)
ln₍d₁₎ = rand(ln_d₁_prior, num_samples)
ln₍s₂₎ = rand(ln_s₂_prior, num_samples)

for i in 1:num_samples
    # p = [ln₍k₆₎[i], ln₍d₁₎[i], ln₍s₂₎[i]] .|> exp 
    p = [ln₍k₆₎[i]] .|> exp

    re_p, u0 = repack_params(updateParams1(p...))
    predictions = solve(problem; p=re_p, saveat=0.1)
    simul_matrix[i,:] = predictions[4,:] + predictions[5,:]
end

my_plot = plot(0:0.1:27.0, simul_matrix'; legend=false, color="#BBBBBB", alpha=0.3)
plot!(
    my_plot,
    0:0.1:27, percentile.(eachcol(simul_matrix), 5);
    fillrange =  percentile.(eachcol(simul_matrix), 95),
    fillalpha = 0.3,
    alpha = 0
    )

plot!(my_plot, 0:0.1:27.0, median.(eachcol(simul_matrix)))

# data = readdlm("Data/fakeOde2/trajectories-0.csv", ',')
# plot!(my_plot, 0:0.1:27, exp.(data[1,:]))
display(my_plot)

    