using DifferentialEquations
using Plots: plot
include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/Differential/ode_restricted.jl")
include("../../Model/treatments_lib.jl")

pred = solve(create_problem(
    max_day=100.0,
    treatment=CBD_IL_12_ver7,
    model="takuya"
    ))
display(plot(pred.t, pred[4,:] + pred[5,:]))

