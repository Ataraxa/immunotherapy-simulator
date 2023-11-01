include("../Library/optimisation_lib.jl")
include("../Model/ode_model.jl")

what_to_plot = "ga_res"

# Paramer vector
if what_to_plot == "ga_res"
    params = load_object("Res/ga_res.jld2")
    params = Evolutionary.minimizer(params)
    global p = params[1:21]
    global u0 = [params[22:end]; 0]
elseif what_to_plot == "benchmark"
    u0, p = get_default_values()
end

fitness([p; u0[1:end-1]])