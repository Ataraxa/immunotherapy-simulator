using Turing 
using DelimitedFiles: readdlm

include("Model/bayesian_hierarchical_model.jl")
include("Model/dde_problem.jl")

# Settings for inference 
time_series_filename = "extensive_cbd7"
is_local_machine = true 

# Create problem object
prob_immune_resp = create_dde_problem()

# Extract data 
data = readdlm()
num_series = sizeof(data)[1]

# Fit data to statistical model 
fitted = fit_hierarchical(data, prob_immune_resp, num_series)

if !is_local_machine 
    chain_dde = Turing.sample(fitted, NUTS(0.7), MCMCDistributed(), 1000, 3; 
        progress=false)

    # Create new filename
    i = 0
    filename = "validation_chain-$i.h5"
    while isfile(filename)
        i+=1
    end
    
    # Save MCMC chain
    h5open("Res/$filename", "w") do f 
        write(f, chain_dde)
    end 
end
