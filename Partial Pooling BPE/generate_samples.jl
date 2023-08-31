# Import standard packages
using Turing
using DifferentialEquations
using StatsPlots: plot, scatter!
using LinearAlgebra
using Random

# Import custom libraries
include("Library/model_lib.jl")
include("Library/stat_lib.jl")

# Settings of the sample generator
Random.seed!(14);
num_traj = 10
do_print = false

# Load default values
u0, p = get_default_values()
tspan = (0.0, 30.0)

for traj in 1:num_traj
    # Overwrite sensitive parameters
    p[2][6] = binorm(0.5, 0.4, 0.8, 0.1, 0.1) # k6 
    p[3][1] = binorm(0.5, 0.4, 0.8, 0.1, 0.1) # d1
    p[4][2] = binorm(0.5, 0.4, 0.8, 0.1, 0.1) # s2 

    h(p, t; idxs::Int) = 1.0 

    prob_dde = DDEProblem(immune_response, u0, h, tspan, p)
    sol_dde = solve(prob_dde; saveat=0.1)
    
    if do_print
        # Plot simulation and noisy observations.
        plot(sol_dde)
        scatter!(sol_dde.t, ddedata[4,:]; color=[1 2], label="")
    end
    
end