"""
Main file where the Bayesian analysis is conducted, according to the proposal.
This file should only be run on the HPC, as it is very heavy.
"""

import Turing
using DifferentialEquations
using DelimitedFiles
using StatsPlots: plot, scatter!
using HDF5
using MCMCChainsStorage

include("Model/bayesian_model.jl")
include("Model/ode_model.jl")

# Settings for inference
traj2fit = "Data/trajectories-average.csv"
is_local_machine = true

# 1 - Create a problem object
u0, p = get_default_values()
t_span = (0.0, 27.0)
h(p, t; idxs::Int) = 0.0
prob_immune_resp = DDEProblem(immune_response, u0, h, t_span, p)

# 2 - Extract data from pre-generated files
data = readdlm(traj2fit, ',')
tracked_tumour_vol = data[5, :] + data[6, :]
approx_sol = Array(tracked_tumour_vol) + 10.0 * randn(size(tracked_tumour_vol)[1])

# seleted_days = 1:trunc(Int, size(tracked_tumour_vol)[1]/15):size(tracked_tumour_vol)[1]
# display(plot(data[1, :], tracked_tumour_vol; color=1))
# scatter!(data[1, seleted_days], tracked_tumour_vol[seleted_days]; color=1)
# scatter!(data[1,seleted_days], approx_sol[seleted_days]; color=2)

# 3 - Fit model to data 
model_dde = fit_immune_resp(approx_sol, prob_immune_resp)

if !is_local_machine
    chain_dde = sample(model_dde, NUTS(0.65), MCMCSerial(), 1000, 3; progress=false)

    h5open("Res/validation_chain.h5", "w") do f 
        write(f, chain_dde)
    end
end

return 0