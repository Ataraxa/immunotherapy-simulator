"""
Script to run a validation protocol on the hierachical Bayesian model
It first creates a matrix of data, and fits it to the desired model.
"""

import Turing
using DotEnv
using TickTock
using HDF5
using MCMCChains
using MCMCChainsStorage

include("../Model/Bayesian/anarchical.jl")
include("../Model/Differential/dde_to_bayesian.jl")

include("../Library/data_extractor.jl")
include("../Library/validation_lib.jl")
include("../Library/optimisation_lib.jl")

tick()
println("Code started running")

# Settings for inference
DotEnv.config() # Loads content from .env file
step_size = parse(Float64, ENV["STEP_SIZE"])
n_iters = (length(ARGS) >= 2) ? parse(Int64, ARGS[1]) : 1000
n_threads = (length(ARGS) >= 2) ? parse(Int64, ARGS[2]) : 1
init_leap = (length(ARGS) >= 3) ? parse(Float64, ARGS[3]) : 0.65
upper1 = (length(ARGS) >= 5) ? parse(Int64, ARGS[4]) : 5
upper2 = (length(ARGS) >= 5) ? parse(Int64, ARGS[5]) : 20

# Create a problem object
prob_immune_resp = restricted_dde_space()

# Extract validation data, select days and add random noise
selected_days = cbd_il_12["active_days"]
num_experiments = 1
data_matrix = cbd_il_12["mean"]

# Fit model to data 
model_dde = fit_anarchical(data_matrix, prob_immune_resp, 
    step_size, selected_days, upper1, upper2)

# This is where the heavy computations come in - so reserved for HPC
pre_sampling = peektimer()
println("Starting sampling ($pre_sampling seconds since last step)")

if ENV["MACHINE_TYPE"] == "hpc"
    println("Going into the remote computing branch")
    chain_dde = Turing.sample(model_dde, NUTS(init_leap), MCMCThreads(), n_iters, n_threads; progress=false)
elseif ENV["MACHINE_TYPE"] == "local" 
    println("Going into the local computing branch")
    chain_dde = Turing.sample(model_dde, NUTS(1.0), MCMCThreads(), 100, 2; progress=false)
end

sampling_time = peektimer() - pre_sampling
println("Finished sampling ($sampling_time seconds since last step)")

# Create new filename
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "$machine-anarchical-$file_i.h5"
while isfile("Res/$filename")
    global file_i+=1
    global filename = "$machine-anarchical-$file_i.h5"
end

# Save MCMC chain
h5open("Res/$filename", "w") do f 
    write(f, chain_dde)
end

# Write in log file
summary = "Summary for $filename: n_iters=$n_iters | n_threads=$n_threads | input_leap=$init_leap | bounds=$upper1,$upper2 \n"
# open("Res/log-$machine.txt", "a") do f 
#     write(f, summary)
# end

# End of scripts log
println("File saved successfully @$filename")
println(summary)