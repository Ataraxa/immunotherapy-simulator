#= FILE DESCRIPTION

Script used to generate fake data purely used for the validation step.
It simply uses a model and a parameter structure to predict the time evolution
of all relevant quantities.

It uses the COMPACT protocol for saving the output: only one file, which 
contains a matrix. Experiments in row direction, 271 time points in column
direction.
=#

using Distributions
using DelimitedFiles

include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/treatments_lib.jl")
include("../../CommonLibrary/struct_manipulation.jl")

### Create problem object
prob = create_problem()

### Generate predictions
predictions = solve(prob; saveat=0.1)
vol = predictions[4,:] + predictions[5,:]

### Add random noise to simulation 
stacked_vol = repeat(vol, 1, 10)'
noise = rand(Normal(0, 0.1), 10, 271)
noisy_vol = log.(stacked_vol) .+ noise

### Save predictions to file
file_i = 0
filename = "trajectories-$file_i.csv"
while isfile("Data/fakeOde2/$filename")
    global file_i+=1
    global filename = "trajectories-$file_i.csv"
end
writedlm("Data/fakeOde2/$filename", noisy_vol, ',')


