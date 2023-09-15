# Import standard packages
using Turing
using DifferentialEquations
using DelimitedFiles
using StatsPlots: plot, scatter!
using LinearAlgebra
using Random
using Distributions
using DotEnv

# Import custom libraries
include("../Model/ode_model.jl")
include("../Model/stat_lib.jl")

# Settings of the sample generator
DotEnv.config()
Random.seed!(14);
num_traj = 10
do_print = true
val_set = ENV["VALIDATION_SET"]

# Load default values
u0, p = get_default_values()
tspan = (0.0, 30.0)

for traj in 1:num_traj

    # Overwrite sensitive parameters
    # p[2][6] = binorm(0.5, 0.2, 0.8, 0.1, 0.1) # k6 
    # p[3][1] = binorm(0.5, 8, 14, 2, 2) # d1
    # p[4][2] = binorm(0.5, 0.1, 0.5, 0.1, 0.1) # s2 

    p[2][6] = rand(truncated(Normal(0.5, 0.3); lower=0)) # k6 
    p[3][1] = rand(truncated(Normal(11, 3); lower=0)) # d1
    p[4][2] = rand(truncated(Normal(0.3, 0.2); lower=0)) # s2 

    h(p, t; idxs::Int) = 1.0 

    prob_dde = DDEProblem(immune_response, u0, h, tspan, p)
    sol_dde = solve(prob_dde; saveat=0.1)

    # Add random (experimental) noise to simulation 
    training_data = Array(sol_dde) + 10 * randn(size(sol_dde))
    
    temp = [sol_dde, training_data]
    if do_print
        # Plot simulation and noisy observations.
        plot(sol_dde, show=true)
        scatter!(sol_dde.t, training_data[4,:]; color=[1 2], label="")
    end
    
    # Save trajectories to .csv file
    matrix_sol = reshape(reduce(hcat, sol_dde), 5, 301)
    writedlm("Data/$val_set/trajectories-$traj.csv", [sol_dde.t'; matrix_sol], ',')
    
    # Save sensitive params to .csv file 
    to_csv = [p[2][6], p[3][1], p[4][2]]
    writedlm("Data/$val_set/params-$traj.csv", to_csv, ',')
end