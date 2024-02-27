#= File Description

This script is used to fit a Bayesian model to a set of manually-
generated data. This is useful to check model identifiability. 
=#
using Dates
using Distributions
using DotEnv
using HDF5
using JLD
using Match
using MCMCChains
using MCMCChainsStorage
using Pipe

include("../../CommonLibrary/data_extractor.jl")
include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")
include("../../Model/Bayesian/individual_restricted.jl")

### Script Settings
DotEnv.config() # Loads content from .env file
step_size = parse(Float64, ENV["STEP_SIZE"])
n_iters         = (length(ARGS) >= 1) ? parse(Int64,   ARGS[1]) : 10
n_threads       = (length(ARGS) >= 2) ? parse(Int64,   ARGS[2]) : 2
num_experiments = (length(ARGS) >= 3) ? parse(Int64,   ARGS[3]) : 1
model           = (length(ARGS) >= 4) ?               (ARGS[4]) : "odeNfullyObs"
prior_distro    = (length(ARGS) >= 5) ?               (ARGS[5]) : "Cauchy"
inform_priors   = (length(ARGS) >= 6) ? parse(Int64,   ARGS[6]) : 0
data_set        = (length(ARGS) >= 7) ? parse(Int64,   ARGS[7]) : 1
prior_acc       = (length(ARGS) >= 8) ? parse(Float64,(ARGS[8])) : 1.0
log_norm        = (length(ARGS) >= 9) ?               (ARGS[9]) : "loga"
σ_err          = (length(ARGS) >= 10) ? parse(Float64,(ARGS[10])) : 1.0

# Manual settings
path = "Data/fake_data"
var_idx = [11,12,21] # For immunotherapy
base = christian_true_params
# var_idx = [1,2,3,4] # For Lotka-Volterra check
# base = [1.5, 1., 3., 1.]

### Settings autoloading
# open("$path/log.txt") do f 
#     lines = readlines(f)
#     for line in lines
#         if string(line[1]) == string(data_set)
#             rhs = split(line, '>')[2]
#             global unparsed_settings = strip.(split(rhs, '|'))
#             break
#         end
#     end
# end
# space = unparsed_settings[1]
# σ_err = 1. # We assume we know the noise level
# log_norm = unparsed_settings[3]
println("σ=$(σ_err) | space=$(space) | transform=$(log_norm)")

### Generate priors 
open("Data/fake_data/params.txt") do f 
    lines = readlines(f)
    for line in lines 
        if string(line[1]) == string(data_set)
            rhs = @pipe split(line, '>')[2]
            rhs2 = strip.(split(rhs, '|'))
            rhs3 = filter(x -> x != "N/A", rhs2) # List of parameters
            global rhs4 = parse.(Float64, rhs3)
            break
        end
    end
end
base[var_idx] .= rhs4
distro = @pipe Symbol(prior_distro) |> getfield(Main, _) # Convert str to distro
priors_vec = gen_priors(distro, prior_acc, Bool(inform_priors); base)
priors_vec = priors_vec[var_idx]

### Main 
# Data Extraction 
# Please refer to convention file to understand the format of csv files
data_mat = load("$path/trajectories-$data_set.jld", "M")
# data_mat = data_mat[:,[1, 20, 30, 40, 50, 60, 70, 80]]
selected_days = [0,7,8,9,11,14,17,20]
data_mat = data_mat[:, selected_days*trunc(Int, 1/step_size) .+ 1,:] #slice pred

# Problem Definition 
problem = create_problem(model=model)

model_args = [data_mat, problem, selected_days, step_size, σ_err,
    base, log_norm, priors_vec, var_idx]
println.(model_args)

fitted_model = fit_individual_restricted3(model_args...; 
    num_experiments = num_experiments)

# fitted_model = fit_lotkaVolterra(model_args...; 
#     num_experiments = num_experiments)

# Sample from Posterior
# println("HERE: $(length(parsed_vec))")
if ENV["MACHINE_TYPE"] == "hpc"
    println("Going into the remote computing branch")
    chain_dde = Turing.sample(fitted_model, NUTS(), MCMCThreads(), n_iters, 
        n_threads; progress=false)
elseif ENV["MACHINE_TYPE"] == "local" 
    println("Going into the local computing branch")
    chain_dde = Turing.sample(fitted_model, NUTS(), MCMCThreads(), n_iters, 
        n_threads; progress=false)
end
display(gelmandiag(chain_dde))


### Enf-of-script log
# Create new filename
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "$machine-individual-$num_experiments-$file_i.h5"
while isfile("Results/$filename")
    global file_i+=1
    global filename = "$machine-individual-$num_experiments-$file_i.h5"
end

# Save MCMC chain
h5open("Results/$filename", "w") do f 
    write(f, chain_dde)
end

# Write in log file
summary = "Summary for $filename: \n 
    MCMC parameters:       n_iters=$n_iters | n_threads=$n_threads \n
    Parameter space:       space=$space  \n
    Priors:                type=$prior_distro | std=$prior_acc | isinform=$inform_priors \n
    Likelihood parameters: model=$model | σ_noise=$σ_err \n 
    Pooling:               n_exp=$num_experiments \n
    Dataset:               set=$data_set \n 
    Transformation:        transform=$log_norm \n \n"

open("Results/log-$machine.txt", "a") do f 
    write(f, summary)
end

println("File saved successfully @$filename ($(now()))")
println(summary)