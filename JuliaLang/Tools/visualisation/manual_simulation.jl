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

plot!(my_plot, 0:0.1:27.0, exp.(data[1:3, :]'); label="Training Data")


savefig(my_plot, "prout.pdf")
