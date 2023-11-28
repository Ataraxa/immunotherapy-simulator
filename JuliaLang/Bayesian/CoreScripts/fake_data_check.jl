#= File Description

This script is used to fit a Bayesian model to a set of manually-
generated data. This is useful to check model identifiability. 
=#
using Pipe
using HDF5
using Dates
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
n_iters         = (length(ARGS) >= 1) ? parse(Int64,   ARGS[1]) : 1000
n_threads       = (length(ARGS) >= 2) ? parse(Int64,   ARGS[2]) : 1
num_experiments = (length(ARGS) >= 3) ? parse(Int64,   ARGS[3]) : 1
model           = (length(ARGS) >= 4) ?               (ARGS[4]) : "takuya"
input_distro    = (length(ARGS) >= 5) ?               (ARGS[5]) : "cauchy"
inform_priors   = (length(ARGS) >= 6) ?               (ARGS[6]) : "true"
data_set        = (length(ARGS) >= 7) ? parse(Int64,   ARGS[7]) : 2

### Settings autoloading
open("Data/fakeData/log.txt") do f 
    lines = readlines(f)
    for line in lines
        if string(line[1]) == string(data_set)
            rhs = split(line, '>')[2]
            global unparsed_settings = strip.(split(rhs, '|'))
            break
        end
    end
end
space = unparsed_settings[2]
σ_likelihood = parse(Float64, unparsed_settings[3])
log_norm = unparsed_settings[4]
println("σ=$(σ_likelihood) | space=$(space) | log_norm=$(log_norm)")

### Generate priors 
open("Data/fakeData/params.txt") do f 
    lines = readlines(f)
    for line in lines 
        if string(line[1]) == string(data_set)
            rhs = split(line, ':')[2]
            rhs2 = strip.(split(rhs, '|'))
            rhs3 = filter(x -> x != "N/A", rhs2)
            # rhs = @pipe split(line, ':')[2] |> strip.(split(_, '|')) |> filter(x -> x != "N/A", _)
            global parsed_vec = [parse(Float64, x) for x in rhs3]
            break
        end
    end
end
println(parsed_vec)

distro = @pipe Symbol(input_distro) |> getfield(Main, _)
ip = parse(Bool, inform_priors)
prior_vec = [distro(ip ? par : 0, 1) for par in parsed_vec]
println(prior_vec)

### Main 
# Data Extraction 
# /!\ Please refer to convention file to understand the format of csv files
selected_days = [0,7,8,9,11,14,17,20]
data_mat = readdlm("Data/fakeData/trajectories-$data_set.csv", ',')
data_mat = data_mat[:, selected_days*trunc(Int, 1/step_size) .+ 1] # slice pred
# scatter(selected_days, data_mat')

# Problem Definition 
problem = create_problem(model=model)

# Fit Model to Data 
@match input_distro begin
    "normal" => global distro = Normal(0,1)
    "cauchy" => global distro = Cauchy(0,1) 
end

model_args = [data_mat, problem, selected_days, step_size, σ_likelihood,
    log_norm, prior_vec]

@match space begin
    "full" => global fitted_model = fit_individual_full(model_args...)
    "restr1" => global fitted_model = fit_individual_restricted1(model_args...; 
        num_experiments = num_experiments)
    "restr3" => global fitted_model = fit_individual_restricted3(model_args...; 
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
    Likelihood parameters: model=$model | σ_noise=$σ_likelihood \n 
    Inference parameters:  n_exp=$num_experiments | distro=$input_distro \n
    Dataset:               set=$data_set \n 
    Transformation:        transform=$log_norm \n \n"

open("Results/log-$machine.txt", "a") do f 
    write(f, summary)
end

println("File saved successfully @$filename ($(now()))")
println(summary)