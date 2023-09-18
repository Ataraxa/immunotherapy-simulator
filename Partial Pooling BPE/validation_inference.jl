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

"Need to wrap whole code within a dummy function for parallel processing"
function multithread_wrapper(create_dde_problem, read_data, fit_unimodal_hierarchical)


    # Settings for inference
    DotEnv.config() # Loads content from .env file
    step_size = parse(Float64, ENV["STEP_SIZE"])
    # println(ENV["MACHINE_TYPE"])

    # Create a problem object
    prob_immune_resp = create_dde_problem()

    # Extract validation data, select days and add random noise
    selected_days = [0,7,8,9,11,14,17,20]
    num_experiments = 1
    data_matrix = read_data(selected_days, num_experiments, step_size)

    # Fit model to data 
    # model_dde = fit_unimodal_hierarchical(data_matrix, prob_immune_resp, num_experiments, 
    #     step_size, selected_days)

    # This is where the heavy computations come in - so reserved for HPC
    if ENV["MACHINE_TYPE"] == "hpc"
        println("Going into the remote computing branch")
        chain_dde = Turing.sample(model_dde, NUTS(0.65), MCMCDistributed(), 1000, 3; progress=false)
    elseif ENV["MACHINE_TYPE"] == "local" 
        println("Going into the local computing branch")
        chain_dde = Turing.sample(fit_unimodal_hierarchical(data_matrix, prob_immune_resp, num_experiments, 
        step_size, selected_days), SMC(), MCMCDistributed(), 10, 1; progress=false)
    end

    # Create new filename
    file_i = 0
    machine = ENV["MACHINE_TYPE"]
    filename = "$machine-validation_chain-$file_i.h5"
    while isfile("Res/$filename")
        file_i+=1
    end

    # Save MCMC chain
    h5open("Res/$filename", "w") do f 
        write(f, chain_dde)
    end

    println("File saved successfully!")
end

multithread_wrapper(create_dde_problem, read_data, fit_unimodal_hierarchical)