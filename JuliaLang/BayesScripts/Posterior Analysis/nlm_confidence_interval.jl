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

chain = h5open("Res/hpc-validation_chain-10.h5", "r") do f
    read(f, Chains)
end

df = nlm_ci(chain[:,:,[1,2,4,5]])
