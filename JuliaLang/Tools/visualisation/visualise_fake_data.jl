using Plots: plot, plot!
using DelimitedFiles: readdlm

include("../../Model/Differential/ode_core.jl")

data = readdlm("Data/fakeOde2/trajectories-1.csv", ',')

# Plot Baseline
base_pb = create_problem()
sol = solve(base_pb; saveat=0.1)
my_plot = plot(0:0.1:27.0, (sol[4,:] + sol[5,:]))

selected_days = [0,7,8,9,11,14,17,20]
plot!(my_plot, selected_days, exp.(data)[:,selected_days*trunc(Int, 1/0.1) .+ 1]')
display(my_plot)