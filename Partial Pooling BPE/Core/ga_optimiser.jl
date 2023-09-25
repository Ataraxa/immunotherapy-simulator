using Evolutionary
using DelimitedFiles
using DifferentialEquations

include("../Model/ode_model.jl")
include("../Library/optimisation_lib.jl")

# Fetch known fit parameters for boxConstraints
u0, p = get_default_values() # u0 has 4 params, p has 21 params
target = [p; u0[1:end-1]] # 1x25 vector
width_factor = 0.5
lower = (1 - width_factor) .* target  # Add unknown params at the end
upper = (1 + width_factor) .* target

# Perform optimisation
params = Evolutionary.optimize(
    fitness,
    BoxConstraints(lower, upper),
    ones(25),
    GA(),
    Evolutionary.Options(parallelization=:thread)
)

# Save result to file
writedlm("Res/ga_results.csv", params, ',')