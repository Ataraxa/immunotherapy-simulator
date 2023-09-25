"""
Script to run a validation protocol on the hierachical Bayesian model
It first creates a matrix of data, and fits it to the desired model.
"""

import Turing
using DotEnv
using HDF5
using MCMCChains
using MCMCChainsStorage

include("../Model/dde_to_bayesian.jl")
include("../Library/data_extractor.jl")
include("../Model/bayesian_unimodal_hierarchical.jl")


# Settings for inference
DotEnv.config() # Loads content from .env file
step_size = parse(Float64, ENV["STEP_SIZE"])

# Create a problem object
prob_immune_resp = restricted_dde_space()

# Extract validation data, select days and add random noise
selected_days = [0,7,8,9,11,14,17,20]
num_experiments = 1
data_matrix = read_data(selected_days, num_experiments, step_size)

# Fit model to data 
model_dde = fit_unimodal_hierarchical(data_matrix, prob_immune_resp, num_experiments, 
    step_size, selected_days)

# This is where the heavy computations come in - so reserved for HPC
if ENV["MACHINE_TYPE"] == "hpc"
    println("Going into the remote computing branch")
    chain_dde = Turing.sample(model_dde, NUTS(0.65), MCMCThreads(), 1000, 3; progress=false)
elseif ENV["MACHINE_TYPE"] == "local" 
    println("Going into the local computing branch")
    chain_dde = Turing.sample(model_dde, NUTS(0.65), 5; progress=false)
end

# Create new filename
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "$machine-validation_chain-$file_i.h5"
while isfile("Res/$filename")
    global file_i+=1
    filename = "$machine-validation_chain-$file_i.h5"
end

# Save MCMC chain
h5open("Res/$filename", "w") do f 
    write(f, chain_dde)
end

println("File saved successfully!")