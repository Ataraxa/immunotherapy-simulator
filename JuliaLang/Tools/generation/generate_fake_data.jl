using Distributions

using JLD
using Pipe
using Match
using Plots: plot, plot!, scatter!

include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")

# Settings
var_params_idx = [11, 12, 21]
noise_std = 1.
noise = "none"
model_name="takuya"

do_overwrite_prev = false

# Sample from prior 
priors = info_p[var_params_idx]
sampled_set = [rand(prior) for prior in priors] .|> exp
param_vector = copy(christian_true_params)
param_vector[var_params_idx] .= sampled_set

# Create and solve 
problem = create_problem(; 
    model=model_name,
    params=param_vector
    )
pred = solve(problem; saveat=0.1) # 5x271 matrix
# v = pred[4,:]+pred[5,:]
v = pred[4,:]
combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
stacked_pred = repeat(combined_pred, 1, 1, 10)

# Add noise
@match noise begin 
    "additive" => begin
        noise = rand(Normal(0, 1)*noise_std, 4, size(stacked_pred)[2], 10)
        global noisy_growth = stacked_pred .+ noise
    end 

    "logn_noise" => begin
        noise = rand(Normal(0, 0.2)*noise_std, 4, size(stacked_pred)[2], 10)   
        global noisy_growth = (log.(stacked_pred) + noise) .|> exp
    end

    "none" => begin
        global noisy_growth = stacked_pred
    end
end

### Rectify data 
function data_rectifier(x, thresh)
    rec_x = (x < thresh) ? thresh : x
    return rec_x
end
data_rectifier.(noisy_growth, 1e-6)

# Plotting result
num_series = 4
series = ["IFNγ", "CD8+", "PD-1", "Tumour volume"]
layout = plot(layout=(2,2))
for i in 1:4
    plot!(layout, pred.t, noisy_growth[i,:,1:num_series];
    sp=i,
    xlabel="Time (day)",
    title=series[i],
    frame=:box)
end
display(layout)

### Save predictions to file
file_i = 0
path = "Data/fakeDataNew"
filename = "trajectories-$file_i.jld"
while isfile("$path/$filename")
    global file_i+=1
    global filename = "trajectories-$file_i.jld"
end
if do_overwrite_prev
    filename = "trajectories-$(file_i-1).jld"
    global file_i -= 1
end
save("$path/$filename", "M", noisy_growth)

summary = "$(file_i) -> $(size(var_params_idx)) | $(noise_std) | $(noise) | $(model_name)\n"
open("$path/log.txt", "a") do f 
    write(f, summary)
end

open("$path/params.txt", "a") do f 
    write(f, "$(file_i) -> "*join(string.(sampled_set), " | ")*"\n")
end

println("Successfully created file at $(filename)")
