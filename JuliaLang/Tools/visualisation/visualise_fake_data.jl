using Plots: plot, plot!, vline!,scatter
using DelimitedFiles: readdlm

include("../../Model/Differential/ode_core.jl")

data = readdlm("Data/fakeData/trajectories-1.csv", ',')

# Plot Baseline
# base_pb = create_problem()
# sol = solve(base_pb; saveat=0.1)
# my_plot = plot(0:0.1:27.0, (sol[4,:] + sol[5,:]))

# selected_days = [0,7,8,9,11,14,17,20]
# my_plot = scatter(selected_days, data[1:5,:]')
# vline!([0,7,8,9,11,14,17,20]; alpha=0.5, legend=false)
plot(0:0.1:27, data')
# display(my_plot)