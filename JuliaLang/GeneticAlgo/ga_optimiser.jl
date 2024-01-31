using Evolutionary
using JLD2
using DifferentialEquations
using DotEnv
using Match

include("./fitness_factory.jl")
include("../Model/Bayesian/priors.jl")

DotEnv.config() # Loads content from .env file
n_iters      = (length(ARGS) >= 1) ? parse(Int64,   ARGS[1]) : 100
n_iter_tol   = (length(ARGS) >= 2) ? parse(Int64,   ARGS[2]) : 10
input_model  = (length(ARGS) >= 3) ?                ARGS[3]  : "takuya"

# Fetch known fit parameters for BoxConstraints
@match input_model begin
    "takuya" => begin
        global p  = christian_true_params[1:21]
        global u0 = [christian_true_params[22:end]; 0]
    end

    "w/feedback" => begin
        global p  = [christian_true_params[1:21]; 1] # last element is n1
        global u0 = [christian_true_params[22:end]; 0]
    end
end
target = [p; u0...] # 1x25 or 1x26 vector
width_factor = 0.25
lower = (1 - width_factor) .* target  # Add unknown params at the end
upper = (1 + width_factor) .* target
fitness = create_fitness(model=input_model)

# Perform optimisation
params = Evolutionary.optimize(
    fitness,
    BoxConstraints(lower, upper),
    GA(
        populationSize=200,
        crossoverRate=0.8,
        mutationRate=0.1
    ),
    Evolutionary.Options(
        parallelization=:thread,
        iterations = n_iters,
        successive_f_tol = n_iter_tol
    )
)

# Save results
file_i = 0
machine = ENV["MACHINE_TYPE"]
filename = "$machine-ga_opt-$file_i.jld2"
while isfile("Results/$filename")
    global file_i+=1
    global filename = "$machine-ga_opt-$file_i.jld2"
end

# Save MCMC chain
save_object("Results/$filename", params)

# Write in log file
# summary = "Summary for $filename: n_iters=$n_iters | n_threads=$n_threads | input_leap=$init_leap | n_exp=$num_experiments \n"
# open("Res/log-$machine.txt", "a") do f 
#     write(f, summary)
# end

# End of scripts log
println("File saved successfully @$filename")
