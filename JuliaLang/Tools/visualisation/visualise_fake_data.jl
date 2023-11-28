using Plots: plot, plot!, vline!, scatter!
using DelimitedFiles: readdlm

include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/Differential/ode_restricted.jl")
include("../../CommonLibrary/struct_manipulation.jl")

data_set = 2
data = readdlm("Data/fakeData/trajectories-$(data_set).csv", ',')

### Visualise raw traces 
# scatter(0:0.1:27, data[]')

### Visualise baseline 
params = exp.([-0.0852590858609503, 1.189102819703836, -0.3274023044544697])
rep, _ = repack_params(updateParams3(params...))
base_pb = create_problem()

sol = solve(base_pb; p=rep, saveat=0.1)
pred = sol[4,:] + sol[5,:]
plot(0:0.1:27, pred)

### Visualise experimental (noisy + sliced)

selected_days = [0,7,8,9,11,14,17,20]
# my_plot = scatter(selected_days, data[1:5,:]')
# vline!([0,7,8,9,11,14,17,20]; alpha=0.5, legend=false)
plot!(selected_days, data[1:5, selected_days*trunc(Int, 1/0.1) .+ 1]'; opacity=0.5)
scatter!(selected_days, data[1:5, selected_days*trunc(Int, 1/0.1) .+ 1]')
# display(my_plot)