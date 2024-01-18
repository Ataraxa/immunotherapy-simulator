using DifferentialEquations
using Plots: plot, plot!, savefig, gr
using DelimitedFiles


include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/Differential/ode_restricted.jl")
include("../../Model/treatments_lib.jl")

gr()
params = christian
# params.d1 = exp(6.1)
# params.k6 = exp(-2)
# params.s2 = exp(-0.5)

pred = solve(create_problem(
    max_day=27.0,
    treatment=CBD_IL_12_ver7,
    model="takuya",
    param_struct=params
    ))
my_plot = plot(pred.t, pred[4,:] + pred[5,:]; label="Posterior Prediction")

# data = readdlm("Data/fakeOde2/trajectories-0.csv", ',')

# selected_days = [0,7,10,12,14,16,18,20,25,27]
# plot!(my_plot, selected_days, (data[1, :])[1:8])
plot!(my_plot, CBD_IL_12_ver7.active_days, CBD_IL_12_ver7.mean;
    markershape=:circle)
display(my_plot)

