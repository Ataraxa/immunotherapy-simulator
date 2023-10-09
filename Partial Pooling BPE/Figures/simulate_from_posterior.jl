"""
This file contains a function to plot simulations of the system using parameters
sampled from the posterior distributions.

This is useful to show that the estimated parameters can correctly reproduce the
experimental results.
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot, plot!, scatter!, savefig 
using StatsBase: percentile, mode

include("../Model/Differential/dde_to_bayesian.jl")

"Plots `n_samples` simulations with parameters sampled from the posterior distributions"
function simulate_from_posterior(chain, n_samples=100)
    # Problem initialisation
    prob_dde = restricted_dde_space()

    # Setup symboles to sample from chain
    k6 = namesingroup(chain, :k6)
    d1 = namesingroup(chain, :d1)
    s2 = namesingroup(chain, :s2)
    
    # Plot settings 
    final_plot = plot(; legend=false, xlim=(0, 27), ylim=(0, 200))

    # Sample and simulate from chain
    simulation_matrix = Array{Float64}(undef, n_samples, 271)
    i = 1
    posterior_samples = sample(chain[[:d1, :s2, :k6]], n_samples; replace=false)
    for p in eachrow(Array(posterior_samples))
        params = [exp(p[3]), exp(p[1]), exp(p[2])]
        
        sol = solve(prob_dde, MethodOfSteps(Tsit5()); p=params, saveat=0.1)
        tumour = sol[4,:] + sol[5,:]
        
        simulation_matrix[i,:] = tumour
        i += 1
    end
    # Plot matrix 
    plot!(0:0.1:27, simulation_matrix' ; alpha=0.1, color="#BBBBBB")
    
    # Plot median and 25-75 quartile range
    plot!(
        0:0.1:27, percentile.(eachcol(simulation_matrix), 25);
        fillrange =  percentile.(eachcol(simulation_matrix), 75),
        fillalpha = 0.3,
        alpha = 0
        )
    plot!(0:0.1:27, median.(eachcol(simulation_matrix)); color=:blue)

    # Plot "true" tumour evolution
    params =  [0.49, 11.31, 0.3434]
    sol_dde = solve(prob_dde; p=params, saveat=0.1)
    tumour = sol_dde[4,:] + sol_dde[5,:]
    scatter!([7,10,12,14,16,18,20,25,27], tumour[[7,10,12,14,16,18,20,25,27]*10])
    savefig(final_plot, "Misc/Images/simul.pdf")
    return final_plot
end

if abspath(PROGRAM_FILE) != @__FILE__
    chain = h5open("Res/hpc-anarchical-2.h5", "r") do f
        read(f, Chains)
    end

    a = simulate_from_posterior(chain, 300)
    display(a)
end