"""
Script to run a validation protocol on the hierachical Bayesian model
It first creates a matrix of data, and fits it to the desired model.
"""

import Turing
using DotEnv
using HDF5
using MCMCChains
using MCMCChainsStorage

include("Model/dde_problem.jl")
include("Tools/data_extractor.jl")
include("Model/unimodal_hierarchical_model.jl")


DotEnv.config() # Loads content from .env file

# Settings for inference
do_plot = true
step_size = parse(Float64, ENV["STEP_SIZE"])
println(ENV["IS_REMOTE"])

# 1 - Create a problem object
prob_immune_resp = create_dde_problem()

# 2 - Extract validation data, select days and add random noise
selected_days = [0,7,8,9,11,14,17,20]
num_experiments = 10
data_matrix = read_data(selected_days, num_experiments, step_size)

# 3 - Fit model to data 
model_dde = fit_dummy_hierarchical(data_matrix, prob_immune_resp, num_experiments, 
    step_size, selected_days)

# This is where the heavy computations come in - so reserved for HPC
if ENV["IS_REMOTE"] == "true"
    println("Going into the remote computing branch")
    chain_dde = Turing.sample(model_dde, NUTS(0.65), MCMCDistributed(), 1000, 3; progress=false)

    # Create new filename
    file_i = 0
    filename = "new_validation_chain-$file_i.h5"
    while isfile("Res/$filename")
        file_i+=1
    end
    
    # Save MCMC chain
    h5open("Res/$filename", "w") do f 
        write(f, chain_dde)
    end

    println("File saved successfully!")
else
    println("Going into the local computing branch")
    chain_dde = Turing.sample(model_dde, SMC(), 3; progress=false)

    # Create new filename
    file_i = 0
    filename = "local_validation_chain-$file_i.h5"
    while isfile("Res/$filename")
        file_i+=1
    end

    # Save MCMC chain
    h5open("Res/$filename", "w") do f 
        write(f, chain_dde)
    end

    println("File saved successfully!")
end

return 0