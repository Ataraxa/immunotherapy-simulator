"""
File where all the analysis on posterior distributions is conducted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot

chain = h5open("Res/hpc-validation_chain-9.h5", "r") do f
    read(f, Chains)
end
a = namesingroup(chain, :k6)
display(plot(chain[:,:,1:2]))

# # Gelman diagnostic: PSRFCI under 1.1 indicates good mixing and convergence
display(gelmandiag(chain)[1:end,:])

# # # Autocorrelation
# display(autocorplot(chain))
