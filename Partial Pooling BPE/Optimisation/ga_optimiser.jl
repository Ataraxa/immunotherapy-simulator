using Evolutionary
using JLD2
using DifferentialEquations

include("../Model/Differential/ode_model.jl")
include("./opt_lib.jl")

# Fetch known fit parameters for boxConstraints
u0, p = get_default_values() # u0 has 4 params, p has 21 params
target = [p; u0[1:end-1]] # 1x25 vector
width_factor = 0.25
lower = (1 - width_factor) .* target  # Add unknown params at the end
upper = (1 + width_factor) .* target

# Block to validate the width factor feature
# for i in eachindex(lower)
#     a = lower[i]
#     b = upper[i]
#     c = target[i]
#     println("$a $c $b")
# end

# Perform optimisation
params = Evolutionary.optimize(
    fitness,
    BoxConstraints(lower, upper),
    lower,
    GA(
        populationSize=200,
        crossoverRate=0.8,
        mutationRate=0.1
    ),
    Evolutionary.Options(
        parallelization=:thread,
        iterations = 10_000,
        successive_f_tol=20_000
    )
)

# Save results
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "$machine-ga_opt-$file_i.h5"
while isfile("Res/$filename")
    global file_i+=1
    global filename = "$machine-ga_opt-$file_i.h5"
end

# Save MCMC chain
h5open("Res/$filename", "w") do f 
    write(f, params)
end

# Write in log file
# summary = "Summary for $filename: n_iters=$n_iters | n_threads=$n_threads | input_leap=$init_leap | n_exp=$num_experiments \n"
# open("Res/log-$machine.txt", "a") do f 
#     write(f, summary)
# end

# End of scripts log
println("File saved successfully @$filename")
