include("./opt_lib.jl")
include("../Model/Differential/ode_model.jl")

u0, p = get_temp_test2()
params = [p; u0]
fitness(params)