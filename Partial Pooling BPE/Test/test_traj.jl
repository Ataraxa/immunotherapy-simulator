#=
File to test if data is correctly extracted before being fitted to the Bayesian
model.
=#

using Plots: plot

include("../Tools/data_extractor.jl")

data_matrix = read_data([0,7,8,9,11,14,17,20], 10, 0.1)
plot(data_matrix') # with plot(), each *column* is a line