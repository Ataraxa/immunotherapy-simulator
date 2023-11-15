using DifferentialEquations
using Plots: plot, plot!, savefig, gr
using DelimitedFiles


include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/Differential/ode_restricted.jl")
include("../../Model/treatments_lib.jl")

gr()
params = christian
params.d1 = exp(6.1)
params.k6 = exp(-2)
params.s2 = exp(-0.5)

pred = solve(create_problem(
    max_day=27.0,
    treatment=CBD_IL_12_ver7,
    model="takuya",
    param_struct=params
    ))
my_plot = plot(pred.t, pred[4,:] + pred[5,:]; label="Posterior Prediction")

data = readdlm("Data/fakeOde2/trajectories-0.csv", ',')

selected_days = [0,7,8,9,11,14,17,20]
plot!(my_plot, selected_days, exp.(data[1, :])[selected_days*trunc(Int, 1/0.1) .+ 1])
display(my_plot)


savefig(my_plot, "prout.pdf")
