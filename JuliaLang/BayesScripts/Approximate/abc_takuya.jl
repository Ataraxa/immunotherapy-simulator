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

# Parse external parameters
param_space     = (length(ARGS) >= 1) ? (ARGS[1]) : "large"

# ABC Settings
println(param_space)
DotEnv.config() # Loads content from .env file
true_params = christian_true_params
priors = gen_priors(Cauchy, 1., false)
if param_space == "large"
    param_indices = collect(1:21)
    println("Inference on large parameter space")
    println(param_indices)
elseif param_space == "medium"
    param_indices = collect(11:21)
    println("Inference on medium parameter space")
    println(param_indices)

else
    param_indices = [11, 12, 21]
    println("Inference on small parameter space")
    println(param_indices)
end

s = 0.1
selected_days = [0,7,8,9,11,14,17,20]
priors = priors[param_indices]
reference_data = load("Data/fakeDataNew/trajectories-4.jld", "M")
reference_data = reference_data[:, selected_days*trunc(Int, 1/s) .+ 1, 1]
output_dir = "Results/abc"
overwrite_last = false

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
    if size(combined_pred, 2) != 271
        println("Skipping sitff combination!")
        println(var_params)
        println("______________________________")
        return zeros(1,32)
    end

    sliced_pred = combined_pred[:,selected_days*trunc(Int, 1/s) .+ 1]
    return sliced_pred
end

# Simulation 
n_particles = 500 # Design choice as well
threshold_schedule = [1000., 750., 500., 250.,175., 100., 75., 50., 25., 17.] # Design choice!
# threshold_schedule = [10., 5., 1., 0.5, 0.3, 0.15, 0.08] # Design choice!
population_colors=["#FFCCD4","#FF667D","#FF2F4E", "#D0001F", "#A20018", 
    "#990017","#800013"]
sim_abcsmc_res = SimulatedABCSMC(
    reference_data,
    simulator,
    priors,
    threshold_schedule,
    n_particles; 

    summary_statistic = "keep_all", # Design choice!
    distance_function = euclidean, # Design choice!
    max_iter=100*n_particles,
    write_progress=true)

# Save result object
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "ABC-$machine-$file_i.jld2"
while isfile("$output_dir/$filename")
    global file_i+=1
    global filename = "ABC-$machine-$file_i.jld2"
end
if overwrite_last
    global filename = "ABC-$machine-$(file_i-1).jld2"
end

save_object("$output_dir/$filename", sim_abcsmc_res)
println("Results saved @ $output_dir/$filename")

# plot(sim_abcsmc_res, population_colors=population_colors)