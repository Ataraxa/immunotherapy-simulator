#= File Description

This script is used to fit a Bayesian model to a set of manually-
generated data. This is useful to check model identifiability. 
=#

using HDF5
using Match
using DotEnv
using MCMCChains
using Distributions
using MCMCChainsStorage
using DelimitedFiles: readdlm

include("../../CommonLibrary/data_extractor.jl")
include("../../Model/Differential/ode_core.jl")
include("../BayesModels/individual_restricted.jl")

### Script Settings
DotEnv.config() # Loads content from .env file
step_size = parse(Float64, ENV["STEP_SIZE"])
n_iters         = (length(ARGS) >= 2) ? parse(Int64,   ARGS[1]) : 1000
n_threads       = (length(ARGS) >= 2) ? parse(Int64,   ARGS[2]) : 1
σ_likelihood    = (length(ARGS) >= 3) ? parse(Float64, ARGS[3]) : 1.0
space           = (length(ARGS) >= 4) ?              (ARGS[4])  : "rest1"
model           = (length(ARGS) >= 5) ?              (ARGS[5])  : "takuya"
num_experiments = (length(ARGS) >= 6) ? parse(Int64,   ARGS[6]) : 1
input_distro    = (length(ARGS) >= 7) ?               (ARGS[7]) : "normal"

### Main 
# Data Extraction
selected_days = [0,7,8,9,11,14,17,20]
data_mat = readdlm("Data/fakeOde2/trajectories-0.csv", ',')

# Problem Definition 
problem = create_problem(model=model)

# Fit Model to Data 
@match input_distro begin
    "normal" => global distro = Normal(0, 1)
    "cauchy" => global distro = Cauchy(0, 1) 
end
model_args = [data_mat, problem, selected_days, step_size, σ_likelihood, distro]

println("Settings for infernce:")
println("Distribution: $distro")
println("Number of experiments: $num_experiments")

@match space begin
    "full" => global fitted_model = fit_individual_full(model_args...)
    "rest1" => global fitted_model = fit_individual_restricted1(model_args...; 
        num_experiments = num_experiments)
    "rest3" => global fitted_model = fit_individual_restricted3(model_args...; 
        num_experiments = num_experiments)
end

# Sample from Posterior
if ENV["MACHINE_TYPE"] == "hpc"
    println("Going into the remote computing branch")
    chain_dde = Turing.sample(fitted_model, NUTS(), MCMCThreads(), n_iters, n_threads; progress=false)
elseif ENV["MACHINE_TYPE"] == "local" 
    println("Going into the local computing branch")
    chain_dde = Turing.sample(fitted_model, NUTS(), MCMCThreads(), 10, 2; progress=false)
end

### Enf-of-script log
# Create new filename
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "$machine-individual-$num_experiments-$file_i.h5"
while isfile("Res/$filename")
    global file_i+=1
    global filename = "$machine-validation_chain-$file_i.h5"
end

# Save MCMC chain
h5open("Results/$filename", "w") do f 
    write(f, chain_dde)
end

# Write in log file
summary = "Summary for $filename: n_iters=$n_iters | n_threads=$n_threads | model=$model | space=$space | n_exp=$num_experiments | distro=$input_distro \n"
open("Results/log-$machine.txt", "a") do f 
    write(f, summary)
end

println("File saved successfully @$filename")
println(summary)