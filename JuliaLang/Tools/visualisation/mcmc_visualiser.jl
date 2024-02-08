"""
File where all the analysis on posterior distributions is conducted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot, savefig

chain = h5open("Results/new_batch/L1reg_n1.h5", "r") do f
# chain = h5open("Results/restr3_n10_log.h5", "r") do f
    read(f, Chains)
end
a = namesingroup(chain, :k6)
plotd = plot(chain[:,:,[2,5]])
# plotd = plot(chain[:,:,:])
display(plotd)
# savefig(plotd,"Misc/Images/batch2/takuya_restr12.png")

# # Gelman diagnostic: PSRFCI under 1.1 indicates good mixing and convergence
display(local df_gelman = gelmandiag(chain[:,:,[2,5]]))

# # # Autocorrelation
# display(autocorplot(chain))
