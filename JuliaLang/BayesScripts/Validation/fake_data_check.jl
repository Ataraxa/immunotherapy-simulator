#= File Description

This script is used to fit a Bayesian model to a set of manually-
generated data. This is useful to check model identifiability. 
=#
using JLD
using Pipe
using HDF5
using Dates
using Match
using DotEnv
using MCMCChains
using Distributions
using MCMCChainsStorage

include("../../CommonLibrary/data_extractor.jl")
include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")
include("../../Model/Bayesian/individual_restricted.jl")

### Script Settings
DotEnv.config() # Loads content from .env file
step_size = parse(Float64, ENV["STEP_SIZE"])
n_iters         = (length(ARGS) >= 1) ? parse(Int64,   ARGS[1]) : 10_000
n_threads       = (length(ARGS) >= 2) ? parse(Int64,   ARGS[2]) : 1
num_experiments = (length(ARGS) >= 3) ? parse(Int64,   ARGS[3]) : 1
model           = (length(ARGS) >= 4) ?               (ARGS[4]) : "odeNnon"
prior_distro    = (length(ARGS) >= 5) ?               (ARGS[5]) : "Normal"
inform_priors   = (length(ARGS) >= 6) ? parse(Int64,   ARGS[6]) : 1
data_set        = (length(ARGS) >= 7) ? parse(Int64,   ARGS[7]) : 0
prior_acc       = (length(ARGS) >= 8) ? parse(Float64,(ARGS[8])) : 1.0
space           = (length(ARGS) >= 9) ?               (ARGS[9]) : "auto"

# Manual settings
path = "Data/fakeDataNew"
space_selection = Dict(
    "restr1" => [1],
    "restr3" => [11, 12, 21],
    "full"   =>  [1:25]
)

### Settings autoloading
open("$path/log.txt") do f 
    lines = readlines(f)
    for line in lines
        if string(line[1]) == string(data_set)
            rhs = split(line, '>')[2]
            global unparsed_settings = strip.(split(rhs, '|'))
            break
        end
    end
end
space = unparsed_settings[1]
σ_err = parse(Float64, unparsed_settings[2]) # We assume we know the noise level
log_norm = unparsed_settings[3]
println("σ=$(σ_err) | space=$(space) | transform=$(log_norm)")

### Generate priors 
open("Data/fakeDataNew/params.txt") do f 
    lines = readlines(f)
    for line in lines 
        if string(line[1]) == string(data_set)
            rhs = @pipe split(line, '>')[2]
            rhs2 = strip.(split(rhs, '|'))
            rhs3 = filter(x -> x != "N/A", rhs2) # List of parameters
            global rhs4 = parse.(Float64, rhs3)
            # rhs = @pipe split(line, ':')[2] |> strip.(split(_, '|')) |> filter(x -> x != "N/A", _)
            break
        end
    end
end
base = christian_true_params
base[[11,12,21]] .= rhs4
distro = @pipe Symbol(prior_distro) |> getfield(Main, _) # Convert str to distro
priors_vec = gen_priors(distro, prior_acc, Bool(inform_priors); base)
# priors_vec = stiff_priors
println.(priors_vec[[11,12,21]])

### Main 
# Data Extraction 
# Please refer to convention file to understand the format of csv files
selected_days = [0,7,8,9,11,14,17,20]
data_mat = load("$path/trajectories-$data_set.jld", "M")
data_mat = data_mat[:, selected_days*trunc(Int, 1/step_size) .+ 1,:] #slice pred

# Problem Definition 
problem = create_problem(model=model)

model_args = [data_mat, problem, selected_days, step_size, σ_err,
    log_norm, priors_vec, [11, 12, 21]]

fitted_model = fit_individual_restricted3(model_args...; 
    num_experiments = num_experiments)

# @match space begin
#     "full" => global fitted_model = fit_individual_full(model_args...)
#     "restr1" => global fitted_model = fit_individual_restricted1(model_args...; 
#         num_experiments = num_experiments)
#     "restr3" => global fitted_model = fit_individual_restricted3(model_args...; 
#         num_experiments = num_experiments)
# end

# Sample from Posterior
# println("HERE: $(length(parsed_vec))")
if ENV["MACHINE_TYPE"] == "hpc"
    println("Going into the remote computing branch")
    chain_dde = Turing.sample(fitted_model, NUTS(), MCMCThreads(), n_iters, 
        n_threads; progress=false)
elseif ENV["MACHINE_TYPE"] == "local" 
    println("Going into the local computing branch")
    chain_dde = Turing.sample(fitted_model, NUTS(), MCMCThreads(), 10, 2; 
        progress=false)
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