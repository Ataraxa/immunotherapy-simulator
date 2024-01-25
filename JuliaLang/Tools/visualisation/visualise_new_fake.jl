using JLD

## Settings
num_series = 4
i = -1
data = load("Data/fakeDataNew/trajectories-$i.jld")

## Plotting 
series = ["IFNγ", "CD8+", "PD-1", "Tumour volume"]
layout = plot(layout=(2,2))
for i in 1:4
    plot!(layout, pred.t, data_mat[i,:,1:num_series];
    sp=i,
    xlabel="Time (day)",
    title=series[i],
    frame=:box)
end
display(layout)