using DifferentialEquations
using Plots: plot, plot!, savefig, gr
using DelimitedFiles


include("../../Model/mechanistic_model.jl")
include("../../Model/treatments_lib.jl")
include("../../Model/Bayesian/priors.jl")

gr()
params = christian_true_params
params[[11, 12, 21]] .= [0.750931634423773, 12.12018890305648, 0.2708004008833685]

pred = solve(create_problem(
    max_day=27.0,
    treatment=CBD_IL_12_ver7,
    model="takuya",
    params=params
    ))
v = pred[4,:]+pred[5,:]
combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
# my_plot = plot(pred.t, pred[4,:] + pred[5,:]; label="Posterior Prediction")
series = ["IFNγ", "CD8+", "PD-1", "Tumour volume"]
layout = plot(layout=(2,2))
for i in 1:4
    plot!(layout, pred.t, combined_pred[i,:];
    sp=i,
    xlabel="Time (day)",
    title=series[i],
    frame=:box)
end
display(layout)
# data = readdlm("Data/fakeOde2/trajectories-0.csv", ',')

# selected_days = [0,7,10,12,14,16,18,20,25,27]
# plot!(my_plot, selected_days, (data[1, :])[1:8])
# plot!(my_plot, CBD_IL_12_ver7.active_days, CBD_IL_12_ver7.mean;
#     markershape=:circle)
# display(my_plot)

