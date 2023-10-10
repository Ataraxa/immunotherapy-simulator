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

n_iters   = (length(ARGS) >= 2) ? parse(Int64,   ARGS[1]) : 1000
n_threads = (length(ARGS) >= 2) ? parse(Int64,   ARGS[2]) : 1
init_leap = (length(ARGS) >= 3) ? parse(Float64, ARGS[3]) : 0.65
std_k6    = (length(ARGS) >= 6) ? parse(Float64, ARGS[4]) : 0.3
std_d1    = (length(ARGS) >= 6) ? parse(Float64, ARGS[5]) : 0.3
std_s2    = (length(ARGS) >= 6) ? parse(Float64, ARGS[6]) : 0.3
max1      = (length(ARGS) >= 8) ? parse(Int64,   ARGS[7]) : 20
max2      = (length(ARGS) >= 8) ? parse(Int64,   ARGS[8]) : 20
exp_err   = (length(ARGS) >= 9) ? parse(Float64, ARGS[9]) : 2

# Create a problem object
prob_immune_resp = restricted_dde_space()

# Extract validation data and logtransform it
selected_days = cbd_il_12["active_days"]
num_experiments = 1
data_matrix = log.(cbd_il_12["mean"])

# Fit model to data 
model_dde = fit_anarchical(data_matrix, prob_immune_resp, 
    step_size, selected_days, std_k6, std_d1, std_s2, max1, max2, exp_err)

# This is where the heavy computations come in - so reserved for HPC
pre_sampling = peektimer()
println("Starting sampling ($pre_sampling seconds since last step)")

if ENV["MACHINE_TYPE"] == "hpc"
    println("Going into the remote computing branch")
    chain_dde = Turing.sample(model_dde, NUTS(), MCMCThreads(), n_iters, n_threads; progress=false)
elseif ENV["MACHINE_TYPE"] == "local" 
    println("Going into the local computing branch")
    chain_dde = Turing.sample(model_dde, NUTS(), MCMCThreads(), 100, 1; progress=false)
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
summary = "Summary for $filename: n_iters=$n_iters | n_threads=$n_threads | input_leap=$init_leap | bounds=$max1,$max2 | std=$std_k6,$std_d1,$std_s2,$exp_err \n"
# open("Res/log-$machine.txt", "a") do f 
#     write(f, summary)
# end

# End of scripts log
println("File saved successfully @$filename")
println(summary)