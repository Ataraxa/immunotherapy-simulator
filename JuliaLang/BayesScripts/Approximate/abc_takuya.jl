using GpABC
using DifferentialEquations

include("../../Model/Differential/ode_core.jl")

# ABC Settings
true_params = [

]
priors = [

]
param_indices = [6, 10, 12]
priors = priors[param_indices]

# Solver Settings ?
dde_problem = create_problem()

# Generator Model 
function simulator(var_params)
    params = copy(true_params)
    params[param_indices] .= var_params

    # Split θ into p and u₀
    p = params[1:21]
    u0 = params[22:end]

    tumour_pred = @pipe solve(dde_problem; p=p, u0=u0) |> +(_[4,:],_[5,:])
    

end

