#= FILE DESCRIPTION

Script used to generate fake data purely used for the validation step.
It simply uses a model and a parameter structure to predict the time evolution
of all relevant quantities.
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
noise = rand(Normal(0, 1), 10, 271)
(vol + noise) .|> log 


### Save predictions to file
matrix_sol = reshape(reduce(hcat, predictions), 5, 301)
writedlm("Data/fakeOde/trajectories-1.csv", [predictions.t'; matrix_sol], ',')

# summary = "Summary for $filename: n_iters=$n_iters | n_threads=$n_threads \n"
# open("Results/log-$machine.txt", "a") do f 
#     write(f, summary)
# end

