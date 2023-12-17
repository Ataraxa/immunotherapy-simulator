"""
File where all the analysis on posterior distributions is conducted
"""

using HDF5
using MCMCChains
using MCMCChainsStorage
using StatsPlots: plot, savefig

chain = h5open("Results/hpc-individual-10￼
Home
Questions
Tags
Users
Unanswered
TEAMS
Stack Overflow for Teams – Start collaborating and sharing organizational knowledge.￼Create a free Team Why Teams?
Consistency of Posterior-2.h5", "r") do f
# chain = h5open("Results/restr3_n10_log.h5", "r") do f
    read(f, Chains)
end
a = namesingroup(chain, :k6)
plotd = plot(chain[:,:,[1,2,3,4,5]])
# plotd = plot(chain[:,:,:])
display(plotd)
savefig(plotd,"Misc/Images/batch2/takuya_restr12.png")

# # Gelman diagnostic: PSRFCI under 1.1 indicates good mixing and convergence
display(local df_gelman = gelmandiag(chain[:,:,[3,4,5]]))

# # # Autocorrelation
# display(autocorplot(chain))
