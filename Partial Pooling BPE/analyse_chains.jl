"""
File where all the analysis on posterior distributions is condicted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot

chain = h5open("Res/new_validation_chain-02.h5", "r") do f
    read(f, Chains)
end

plot(chain)