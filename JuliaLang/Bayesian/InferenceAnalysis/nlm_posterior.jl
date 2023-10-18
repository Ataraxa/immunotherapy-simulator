using Plots: plot, plot!,density
using HDF5
using MCMCChains
using MCMCChainsStorage
using Distributions
using StatsPlots
using LaTeXStrings

"Draws posterior distributions for fitted non-multilevel (nml) model"
function draw_nml_posterior(nml_chains)
    plot(nml_chains)
end

chain = h5open("Res/hpc-anarchical-2.h5", "r") do f
    read(f, Chains)
end

draw_nml_posterior(chain)

