include("../Library/optimisation_lib.jl")
include("../Model/ode_model.jl")

u0, p = get_default_values()
fitness([p; u0[1:end-1]])