using Evolutionary
using JLD2
using DifferentialEquations

include("../Model/ode_model.jl")
include("../Library/optimisation_lib.jl")

# Fetch known fit parameters for boxConstraints
u0, p = get_default_values() # u0 has 4 params, p has 21 params
target = [p; u0[1:end-1]] # 1x25 vector
width_factor = 0.1
lower = (1 - width_factor) .* target  # Add unknown params at the end
upper = (1 + width_factor) .* target

# Perform optimisation
params = Evolutionary.optimize(
    fitness,
    BoxConstraints(lower, upper),
    ones(25),
    GA(
        populationSize=200,
        crossoverRate=0.8,
        mutationRate=0.1
    ),
    Evolutionary.Options(
        parallelization=:thread,
        iterations = 10000
        successive_f_tol=2000
    )
)

# Save result to file
save_object("Res/ga_res.jld2", params)