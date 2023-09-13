"""
Script to run a validation protocol on the hierachical Bayesian model
It first creates a matrix of data, and fits it to the desired model.
"""

import Turing
using DifferentialEquations
using DelimitedFiles
using StatsPlots: plot, scatter!
using HDF5
using MCMCChains
using MCMCChainsStorage

include("Model/bayesian_base_model.jl")
include("Model/dde_problem.jl")

# Settings for inference
do_plot = true

# 1 - Create a problem object
prob_immune_resp = create_dde_problem()

# 2 - Extract validation data, select days and add random noise
data_matrix
data = readdlm(traj2fit, ',')
exact_vol = data[5, :] + data[6, :]
approx_vol = Array(exact_vol) + 10.0 * randn(size(exact_vol)[1])
# selected_days = 1:trunc(Int, size(exact_vol)[1]/15):size(exact_vol)[1]
selected_days = [0, 7, 8, 9, 11, 14, 17, 20]

if do_plot
    display(plot(data[1, :], exact_vol; color=1))
    display(scatter!(data[1,selected_days*10 .+ 1], approx_vol[selected_days*10 .+ 1]; color=1))
end 

validation_data = approx_vol[selected_days*10 .+ 1]

# 3 - Fit model to data 
model_dde = fit_immune_resp(approx_vol, prob_immune_resp)

# This is where the heavy computations come in - so reserved for HPC
if !is_local_machine
    chain_dde = Turing.sample(model_dde, NUTS(0.65), MCMCDistributed(), 1000, 2; progress=false)
    
    # Create new filename
    i = 0
    filename = "validation_chain-$i.h5"
    while isfile(filename)
        i+=1
    end
    
    # Save MCMC chain
    h5open("Res/$filename", "w") do f 
        write(f, chain_dde)
    end
end

return 0