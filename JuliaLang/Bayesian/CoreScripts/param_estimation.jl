#= FILE DESCRIPTION

Script used to estimate parameters of the DDE model by fitting it to a time 
series. This is particulary useful to verify that a Bayesian analysis can yield 
the correct results when we already know the correct parameters beforehand. 
=#

import Turing # to run bayesian inference 
import DifferentialEquations
import MCMCChains
import MCMCChainsStorage
import HDF5: h5open
import Match

include("../../Model/Differential/ode_core.jl")
include("../../CommonLibrary/data_extractor.jl")
include("../BayesModels/individual_full_model.jl")
include("../BayesModels/individual_restricted.jl")

#### Settings for inference
# space = either full (21 parameters) or restricted (2-5 params)
# model = takuya, feedbacked, etc...
####

DotEnv.config() # Loads content from .env file
step_size = parse(Float64, ENV["STEP_SIZE"])

n_iters         = (length(ARGS) >= 2) ? parse(Int64,   ARGS[1]) : 1000
n_threads       = (length(ARGS) >= 2) ? parse(Int64,   ARGS[2]) : 1
σ_likelihood    = (length(ARGS) >= 3) ? parse(Float64, ARGS[3]) : 2.0
space           = (length(ARGS) >= 4) ?                ARGS[4]  : "rest1"
model           = (length(ARGS) >= 5) ?                ARGS[5]  : "takuya"

### Setting up the inference 
# DDE Problem
problem = create_problem(model=model)

# Data Extraction 
selected_days = [0,7,8,9,11,14,17,20]
temp = read_data(selected_days, 1, step_size, "fakeOde")
println(temp)
data_vector = log.(temp)

### Run inference
model_args = [data_vector, problem, selected_days, step_size, σ_likelihood]
@match space begin
    "full" => global fitted_model = fit_individual_full(model_args...)
    "rest1" => global fitted_model = fit_individual_restricted1(model_args...)
    "rest3" => global fitted_model = fit_individual_restricted3(model_args...)
end

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
summary = "Summary for $filename: n_iters=$n_iters | n_threads=$n_threads | model=$model | space=$space \n"
open("Results/log-$machine.txt", "a") do f 
    write(f, summary)
end

println("File saved successfully @$filename")
println(summary)