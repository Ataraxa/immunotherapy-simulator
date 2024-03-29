#= File Description

This file contains a script to run a prior predictive check on a model.
User must specify both the priors and likelihood.
=#

using Plots: plot, plot!, savefig, pyplot
# using Turing
using StatsBase: percentile, mode, median
using Distributions
using DelimitedFiles
using Plots.PlotMeasures

include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")


### Script Settings
# pyplot()
num_samples = 1000
simul_matrix = Matrix{Float64}(undef, num_samples, 271) # (undef, row, col)
problem = create_problem()

### Priors 
ln_k₆_prior = truncated(Cauchy(0, 1); lower=-100, upper=4)
ln_d₁_prior = truncated(Cauchy(0, 1); lower=-100,  upper=4)
ln_s₂_prior = truncated(Cauchy(0, 1); lower=-100, upper=4)

### Main
ln₍k₆₎ = rand(ln_k₆_prior, num_samples)
ln₍d₁₎ = rand(ln_d₁_prior, num_samples)
ln₍s₂₎ = rand(ln_s₂_prior, num_samples)

for i in 1:num_samples
    p = [ln₍k₆₎[i], ln₍d₁₎[i], ln₍s₂₎[i]] .|> exp 
    # p = [ln₍k₆₎[i]] .|> exp

    params = christian_true_params
    params[[11, 12, 21]] .= p
    predictions = solve(problem; p=params, saveat=0.1)
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

plot!(my_plot, 0:0.1:27.0, median.(eachcol(simul_matrix)); 
    yaxis=:linear,
    xlabel="Time (days)",
    ylabel="Tumour volume (mm³)",
    colour=:red,
    left_margin=10mm,
    dpi=600)

# data = readdlm("Data/fakeOde2/trajectories-0.csv", ',')
# plot!(my_plot, 0:0.1:27, exp.(data[1,:]))
display(my_plot)

savefig(my_plot, "LOGprout.png")