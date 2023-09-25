"""
File where all the analysis on posterior distributions is condicted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot

chain = h5open("Res/local-validation_chain-0.h5", "r") do f
    read(f, Chains)
end

plot(chain)