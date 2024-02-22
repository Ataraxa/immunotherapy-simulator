using Pipe
using JLD
using JLD2
using Plots
using GpABC
using DotEnv
using Distances: euclidean, sqeuclidean, peuclidean
using Distributions
using DifferentialEquations

include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")
s=0.1
selected_days = [0,7,8,9,11,14,17,20]
reference_data = load("Data/fakeDataNew/trajectories-0.jld", "M")
reference_data = reference_data[:, selected_days*trunc(Int, 1/s) .+ 1, 1]
reference_data = reshape(reference_data, (1, 32))

sliced_pred = zeros(1,32)
# params = copy(christian_true_params)
# params[[11, 12, 21]] .= [1.4104497537492389, 0.7146079984999922, -0.2983212967417186]
# prob = create_problem(; model="takuya")
# pred = solve(prob, p=params, saveat=s)
# v = pred[4,:]+pred[5,:]
# combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
# sliced_pred = combined_pred[:,selected_days*trunc(Int, 1/s) .+ 1]
# sliced_pred = reshape(sliced_pred, (1, 32))

distance = sqrt(sum((sliced_pred - reference_data).^2))



