using JLD
using Plots

include("../../Model/mechanistic_model.jl")
# Settings
num_series = 3
set_index = 7
data_mat = load("Data/fake_data/trajectories-$set_index.jld", "M")
# data_mat = data_mat[:,[1, 20, 30, 40, 50, 60, 70, 80]]

# LV_model
sol = solve(create_problem(model="takuya"))
# Plotting 
series = ["IFNγ", "CD8+", "PD-1", "Tumour volume"]
layout = plot(layout=(2,2))
for i in 1:4
    plot!(layout, 0:0.1:27, data_mat[i,:,1:num_series];
    sp=i,
    xlabel="Time (day)",
    title=series[i],
    frame=:box)
    # plot!(sol.t, sol[i])
end
display(layout)

# Plot simulation and noisy observationsf for Lotka-Volterra
# plot(sol; alpha=0.3)
# plt1 = plot()
# plot!(sol; label=["x" "y"])
# scatter!([0, 2, 3, 4, 5, 6, 7, 8], data_mat'; color=[1 2], label=["x" "y"])
# xlabel!("Time")
# ylabel!("Population Density")
