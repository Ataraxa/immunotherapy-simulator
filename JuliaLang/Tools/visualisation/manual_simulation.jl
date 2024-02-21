using OrdinaryDiffEq
using Plots: plot, plot!, savefig, gr
using DelimitedFiles


include("../../Model/mechanistic_model.jl")
include("../../Model/treatments_lib.jl")
include("../../Model/Bayesian/priors.jl")
println("Starting")
gr()
params1 = christian_true_params
# params1[[11, 12, 21]] .= [0.4037028186937715, 16.470910037502286, 0.6457205405500415]
params1[[11, 12, 21]] .= ([0.6290015460402033, 1.545447608923281, 0.4329922695466099])

problem = create_problem(; model="odeNfullyObs", max_day=27.0)

pred = solve(problem, KenCarp5(); p=params1)
v = pred[4,:]
combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
# pred.t
my_plot = plot(pred.t, pred[4,:] + pred[5,:]; label="Posterior Prediction")
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

