"""
File to debug the Bayesian model
"""

using Turing
using DifferentialEquations
using StatsPlots: plot
using HDF5
using MCMCChainsStorage

include("../Model/bayesian_model.jl")
include("../Model/ode_model.jl")

u0, p = get_default_values()
t_span = (0.0, 30.0)
h(p, t; idxs::Int) = 0.0
prob_immune_resp = DDEProblem(immune_response, u0, h, t_span, p)

test_model = fit_immune_resp(missing, prob_immune_resp)
c = sample(test_model, SMC(), 1000)

h5open("Res/test_chain.h5", "w") do f 
    write(f, c)
end