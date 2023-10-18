using HDF5
using MCMCChains
using MCMCChainsStorage

# Create new filename
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "$machine-validation_chain-$file_i.h5"
while isfile("Res/$filename")
    global file_i+=1
    filename = "$machine-validation_chain-$file_i.h5"
end

# Save MCMC chain
h5open("Res/$filename", "w") do f 
    write(f, chain_dde)
end

println("File saved successfully!")