"""
File where all the analysis on posterior distributions is conducted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot

chain = h5open("Res/local-validation_chain-3.h5", "r") do f
    read(f, Chains)
end
plot(chain)
# Gelman diagnostic: PSRFCI under 1.1 indicates good mixing and convergence
# display(gelmandiag(chain)[1:end,:])

# # Autocorrelation and independent sampling analysis
# autocorplot(chain)
