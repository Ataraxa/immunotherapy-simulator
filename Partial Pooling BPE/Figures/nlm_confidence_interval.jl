using Plots: plot, plot!,density
using HDF5
using MCMCChains
using MCMCChainsStorage
using Distributions
using StatsPlots
using LaTeXStrings

function nlm_ci(chain)
    hpd(chain)
end

chain = h5open("Res/hpc-anarchical-2.h5", "r") do f
    read(f, Chains)
end

nlm_ci(chain)
