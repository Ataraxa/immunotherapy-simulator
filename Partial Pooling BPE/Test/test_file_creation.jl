#=
File that creates dummy empty files. Useful to test functions that save
MCMC chains in .h5 files.
=#
using HDF5
using MCMCChains
using MCMCChainsStorage

println("Starting sampling process")
dummy_chain = Chains(randn(50, 2, 4), [:a, :b])
println("Finished")

file_i = 0
machine = "dummy"
filename = "$machine-validation_chain-$file_i.h5"
while isfile("Res/$filename")
    global file_i+=1
    filename = "$machine-validation_chain-$file_i.h5"
    println("one up!")
end

# Save MCMC chain
h5open("Res/$filename", "w") do f 
    write(f, dummy_chain)
end