using DelimitedFiles

include("../Model/Differential/ode_core.jl")
include("../Model/Differential/ode_params.jl")
include("../Model/treatments_lib.jl")
include("../CommonLibrary/struct_manipulation.jl")

### Create problem object
p, u0 = struct_split(christian)
p = [p; CBD_IL_12_ver7]
u0 = [u0; 0]
h(p, t; idxs::Int) = 0.0
t_span = (0.0, 30.0)
prob = DDEProblem(full_immune_response, u0, h, t_span, p)

### Generate predictions
predictions = solve(prob; p=p, saveat=0.1)

### Add random noise to simulation 
# training_data = Array(predictions) + 10 * randn(size(predictions))

### Save predictions to file
matrix_sol = reshape(reduce(hcat, predictions), 5, 301)
writedlm("Data/fakeOde/trajectories-1.csv", [predictions.t'; matrix_sol], ',')

