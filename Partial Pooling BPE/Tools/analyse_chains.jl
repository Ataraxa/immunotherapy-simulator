"""
File where all the analysis on posterior distributions is conducted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot

chain = h5open("Res/hpc-validation_chain-10.h5", "r") do f
    read(f, Chains)
end
a = namesingroup(chain, :k6)
display(plot(chain[:,:,[1,2,4,5]]))

# # Gelman diagnostic: PSRFCI under 1.1 indicates good mixing and convergence
df_gelman = gelmandiag(chain[:,:,[1,2,4,5]])[1:end,:]

# # # Autocorrelation
# display(autocorplot(chain))
