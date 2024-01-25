using JLD
using Plots: plot, plot!
# Settings
num_series = 4
i = 0
data_mat = load("Data/fakeDataNew/trajectories-$i.jld", "M")

# Plotting 
series = ["IFNγ", "CD8+", "PD-1", "Tumour volume"]
layout = plot(layout=(2,2))
for i in 1:4
    plot!(layout, 0:0.1:27, data_mat[i,:,1:num_series];
    sp=i,
    xlabel="Time (day)",
    title=series[i],
    frame=:box)
end
display(layout)