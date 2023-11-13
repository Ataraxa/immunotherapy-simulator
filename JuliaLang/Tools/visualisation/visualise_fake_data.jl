using Plots: plot, plot!
using DelimitedFiles: readdlm

include("../../Model/Differential/ode_core.jl")

data = readdlm("Data/fakeOde2/trajectories-0.csv", ',')

# Plot Baseline
base_pb = create_problem()
sol = solve(base_pb; saveat=0.1)
my_plot = plot(0:0.1:27.0, (sol[4,:] + sol[5,:]))


plot!(my_plot, 0:0.1:27.0, exp.(data[1, :]))
display(my_plot)