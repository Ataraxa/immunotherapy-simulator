using HDF5
using MCMCChains
using MCMCChainsStorage

# Construct a chain and write it out...
chain = Chains(randn(500, 2, 4), [:a, :b])
h5open("an_hdf5_file.h5", 'w') do f
  write(f, chain)
end

# # ...and we can get it back
# chain = h5open("an_hdf5_file.h5", "r") do f
#   read(f, Chains)
# endusing HDF5
using MCMCChains
using MCMCChainsStorage

# Construct a chain and write it out...
chain = Chains(randn(500, 2, 4), [:a, :b])
h5open("an_hdf5_file.h5", "w") do f
  write(f, chain)
end