
using JLD
using Plots
using GpABC
using Distances
using Distributions
using DifferentialEquations

include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")

# ABC Settings
true_params = christian_true_params
priors = gen_priors(Normal, 1., true)
param_indices = [11, 12, 21]
s = 0.1
selected_days = [0,7,8,9,11,14,17,20]
priors = priors[param_indices]
reference_data = load("Data/fakeDataNew/trajectories-0.jld", "M")
reference_data = reference_data[:, selected_days*trunc(Int, 1/s) .+ 1, 1]

# Solver Settings ?
dde_problem = create_problem(model="takuya")

# Generator Model 
function simulator(var_params)
    var_params = exp.(var_params)
    # println(var_params)
    params = copy(true_params)
    params[param_indices] .= var_params

    # Split θ into p and u₀
    p = params[1:21]
    u0 = [params[22:end]; 0]

    pred = solve(dde_problem; p=p, u0=u0, verbose=true, saveat=0.1)
    v = pred[4,:]+pred[5,:]
    combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))

    # Debug zone
    # println(var_params)
    # println(size(combined_pred))
    # println(size(pred))
    # println("____________________________________")

    sliced_pred = combined_pred[:,selected_days*trunc(Int, 1/s) .+ 1]

    return sliced_pred
end

# Distance functions
mse(x,y) = mean((x-y).^2)

# Simulation 
n_particles = 200
threshold_schedule = [1000., 500., 250., 100., 50., 40., 10.] # Design choice!
population_colors=["#FF2F4E", "#D0001F", "#A20018", "#990017", "#800013"]
sim_abcsmc_res = SimulatedABCSMC(reference_data,
    simulator,
    priors,
    threshold_schedule,
    n_particles; 
    distance_function = euclidean,
    summary_statistic = "keep_all",
    write_progress=false)
plot(sim_abcsmc_res, population_colors=population_colors)