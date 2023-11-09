# Import standard packages
using Turing
using DifferentialEquations
using DelimitedFiles
using StatsPlots: plot, scatter!
using LinearAlgebra
using Random

# Import custom libraries
include("../Model/ode_model.jl")
include("../Model/stat_lib.jl")

# Settings of the sample generator
Random.seed!(14);
num_traj = 10
do_print = true
step_size = 0.1

# Load default values
u0, p = get_default_values()
tspan = (0.0, 27.0)

h(p, t; idxs::Int) = 0.0 

prob_dde = DDEProblem(immune_response, u0, h, tspan, p)
sol_dde = solve(prob_dde; saveat=step_size)

# Add random (experimental) noise to simulation 
training_data = Array(sol_dde) + 10 * randn(size(sol_dde))

temp = [sol_dde, training_data]


# Save trajectories to .csv file
dim0 = size(sol_dde)
matrix_sol = reshape(reduce(hcat, sol_dde), dim0[1], dim0[2])
writedlm("Data/trajectories-average.csv", [sol_dde.t'; matrix_sol], ',')

# Save sensitive params to .csv file 
to_csv = [p[2][6], p[3][1], p[4][2], step_size]
writedlm("Data/params-average.csv", to_csv, ',')


0