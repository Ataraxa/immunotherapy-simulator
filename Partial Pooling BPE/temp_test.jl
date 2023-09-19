a, b, c,
d, e, f = [1, 2,3, 4,5, 6]


# using DifferentialEquations
# using Distributions
# include("Model/dde_problem.jl")

# @model function testmodel()
    
#     problem = create_dde_problem()
#     p = (0.3, 11.0, 0.3)
#     predictions = solve(problem; p=p, saveat=0.1)
#     pred_vol = predictions[4,:] + predictions[5,:]
#     sliced_pred = pred_vol[[0, 7, 9, 11, 14, 21]*trunc(Int, 1/0.1) .+ 1]
#     x ~ Normal(sliced_pred[1], 1)
# end

# sample(testmodel(), SMC(), 10)


# using HDF5
# using MCMCChains
# using MCMCChainsStorage
# using StatsPlots: plot

# chain = h5open("Res/yet_a_test.h5", "r") do f
#     read(f, Chains)
# end

