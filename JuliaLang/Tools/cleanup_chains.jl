"""
File where all the analysis on posterior distributions is conducted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot

filename ="hpc-validation_chain-10"
chain = h5open("Res/$filename.h5", "r") do f
    read(f, Chains)
end

corrected_chains = chain[:,:,[1,2,4,5]]

h5open("Res/$filename-clup.h5", "w") do f 
    write(f, corrected_chains)
end

