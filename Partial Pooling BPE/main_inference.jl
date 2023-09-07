"""
Main file where the Bayesian analysis is conducted, according to the proposal.
This file should only be run on the HPC, as it is very heavy.
"""

import Turing: sample
using DifferentialEquations
using DelimitedFiles
using StatsPlots
using HDF5
using MCMCChainsStorage

include("Model/bayesian_model.jl")
include("Model/ode_model.jl")

# Settings for inference
traj2fit = "Data/trajectories-1.csv"

# 1 - Create a problem object
u0, p = get_default_values()
t_span = (0.0, 30.0)
h(p, t; idxs::Int) = 0.0
prob_immune_resp = DDEProblem(immune_response, u0, h, t_span, p)

# 2 - Extract data from pre-generated files
global data = readdlm(traj2fit, ',')
tracked_tumour_vol = data[5, :] + data[6, :]

# 3 - Fit model to data 
# println(size(tracked_tumour_vol))
model_dde = fit_immune_resp(tracked_tumour_vol, prob_immune_resp)

chain_dde = sample(model_dde, SMC(), MCMCSerial(), 3000, 3; progress=false)

h5open("Res/validation_chain.h5", 'w') do f 
    write(f, chain_dde)
end
