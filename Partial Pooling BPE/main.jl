"""
Main file where the Bayesian analysis is conducted, according to the proposal
"""

import Turing: sample
import DifferentialEquations
import DelimitedFiles

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
tracked_tumour_vol = data[5, :], data[6, :]

# 3 - Fit model to data 
model_dde = fit_immune_resp(tracked_tumour_vol, prob_immune_resp)

chain_dde = sample(model_dde, NUTS(0.65), MCMCSerial(), 300, 3; progress=false)

