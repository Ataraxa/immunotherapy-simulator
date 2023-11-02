"""
Script to estimate the 25 parameters of a tumour DDE model through a Bayesian 
statistical model.
"""

import Turing # to run bayesian inference 
import DifferentialEquations
import MCMCChains
import MCMCChainsStorage
import HDF5: h5open

include("../../Model/Differential/ode_core.jl")
include("../../CommonLibrary/data_extractor.jl")
include("../BayesModels/individual_full_model.jl")

### Settings for inference
DotEnv.config() # Loads content from .env file
step_size = parse(Float64, ENV["STEP_SIZE"])

n_iters         = (length(ARGS) >= 2) ? parse(Int64,   ARGS[1]) : 1000
n_threads       = (length(ARGS) >= 2) ? parse(Int64,   ARGS[2]) : 1
σ_likelihood    = (length(ARGS) >= 3) ? parse(Float64, ARGS[3]) : 2.0

### Setting up the inference 
# DDE Problem
problem = problem_factory()

# Data Extraction 
selected_days = [0,7,8,9,11,14,17,20]
data_vector = log.(read_data(selected_days, 1, step_size, "fakeOde"))

### Run inference 
fitted_model = fit_individual(data_vector, problem, 
    selected_days, step_size, σ_likelihood)

### Sample from Posterior
if ENV["MACHINE_TYPE"] == "hpc"
    println("Going into the remote computing branch")
    chain_dde = Turing.sample(fitted_model, NUTS(), MCMCThreads(), n_iters, n_threads; progress=false)
elseif ENV["MACHINE_TYPE"] == "local" 
    println("Going into the local computing branch")
    chain_dde = Turing.sample(fitted_model, NUTS(), MCMCThreads(), 10, 2; progress=false)
end

### End of script logs 
# Create new filename
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "$machine-individual-$file_i.h5"
while isfile("Res/$filename")
    global file_i+=1
    global filename = "$machine-validation_chain-$file_i.h5"
end

# Save MCMC chain
h5open("Results/$filename", "w") do f 
    write(f, chain_dde)
end

# Write in log file
summary = "Summary for $filename: n_iters=$n_iters | n_threads=$n_threads \n"
open("Results/log-$machine.txt", "a") do f 
    write(f, summary)
end

println("File saved successfully @$filename")
println(summary)