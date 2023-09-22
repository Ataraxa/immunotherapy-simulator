#=
File to run a series of test on hierarchical statistical models
=#
using StatsPlots: plot
include("../Model/dde_problem.jl")
include("../Model/unimodal_hierarchical_model.jl")
include("../Model/bayesian_hierarchical_model.jl")

test_model = fit_unimodal_hierarchical(missing, create_dde_problem(),10,0.1,[0,7,8,9,11,14,17,20])
c = sample(test_model, SMC(), 10)
plot(c)