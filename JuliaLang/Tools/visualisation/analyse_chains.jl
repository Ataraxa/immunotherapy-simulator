"""
File where all the analysis on posterior distributions is conducted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot, savefig

chain = h5open("Results/hpc-individual-0.h5", "r") do f
    read(f, Chains)
end
a = namesingroup(chain, :k6)
plotd = plot(chain)
savefig(plotd,"Misc/Images/batch2/takuya_restr1.png")

# # Gelman diagnostic: PSRFCI under 1.1 indicates good mixing and convergence
# df_gelman = gelmandiag(chain[:,:,:])

# # # Autocorrelation
# display(autocorplot(chain))
