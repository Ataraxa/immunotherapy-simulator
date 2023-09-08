"""
Main file where the Bayesian analysis is conducted, according to the proposal.
This file should only be run on the HPC, as it is very heavy.
"""

import Turing
using DifferentialEquations
using DelimitedFiles
using StatsPlots: plot, scatter!
using HDF5
using MCMCChains
using MCMCChainsStorage

include("Model/bayesian_model.jl")
include("Model/dde_problem.jl")

# Settings for inference
traj2fit = "Data/trajectories-average.csv"
is_local_machine = false

# 1 - Create a problem object
prob_immune_resp = create_dde_problem()

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