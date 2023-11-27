using Distributions
using DelimitedFiles
using Plots: plot, plot!, scatter!
using Match

include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_restricted.jl")

# Settings
do_overwrite_prev = false
distro_name = "cauchy" 
space = "restr1"
noise = "logn_noise"

# Process Settings
@match distro_name begin 
    "cauchy" => global distro = Cauchy(0,1)
    "normal" => global distro = Normal(0,1)
end

# Priors to sample from 
ln_k₆_prior = truncated(distro; lower=-100, upper=0)
ln_d₁_prior = truncated(distro; lower=0, upper=7)
ln_s₂_prior = truncated(distro; lower=-100, upper=0)
σ_prior     = Normal()

# Sample from priors 
ln₍k₆₎ = rand(ln_k₆_prior)
ln₍d₁₎ = rand(ln_d₁_prior)
ln₍s₂₎ = rand(ln_s₂_prior)
σ = 0.3

# Create and solve 
problem = create_problem()

@match space begin 
    "restr1" => begin  
        # p = [ln₍k₆₎] .|> exp
        p = [exp(-0.8733718288575096)]
        global re_p, _ = repack_params(updateParams1(p...))
    end

    "restr3" => begin
        p = [ln₍k₆₎, ln₍d₁₎, ln₍s₂₎] .|> exp

        global re_p, _ = repack_params(updateParams3(p...))
    end
end

# p = [ln₍k₆₎[i], ln₍d₁₎[i], ln₍s₂₎[i]] .|> exp 
# p = [ln₍k₆₎] .|> exp

# re_p, u0 = repack_params(updateParams1(p...))
predictions = solve(problem; p=re_p, saveat=0.1)
mean_growth = predictions[4,:] + predictions[5,:] 
stacked_growth = repeat(mean_growth, 1, 10)'

@match noise begin 
    "norm_noise" => begin
        noise = rand(Normal(0, 1), 10, size(stacked_growth)[2])
        global noisy_growth = stacked_growth .+ noise
    end 

    "logn_noise" => begin
        noise = rand(Normal(0, 0.2), 10, size(stacked_growth)[2])   
        global noisy_growth = (log.(stacked_growth) + noise) .|> exp
    end
end

### Rectify data 
function data_rectifier(x, thresh)
    rec_x = (x < thresh) ? thresh : x
    return rec_x
end
data_rectifier.(noisy_growth, 1e-6)

my_plot = plot(predictions.t, mean_growth)
plot!(my_plot, predictions.t, noisy_growth')
display(my_plot)

### Save predictions to file
file_i = 0
filename = "trajectories-$file_i.csv"
while isfile("Data/fakeData/$filename")
    global file_i+=1
    global filename = "trajectories-$file_i.csv"
end
if do_overwrite_prev
    filename = "trajectories-$(file_i-1).csv"
    global file_i -= 1
end
writedlm("Data/fakeData/$filename", noisy_growth, ',')

summary = "$(file_i) -> $(distro_name) | $(space) | $(σ) | $(noise)\n"
open("Data/fakeData/log.txt", "a") do f 
    write(f, summary)
end

@match space begin 
    "restr1" => global param_list = "$file_i: $ln₍k₆₎ | N/A | N/A\n"
    "restr3" => global param_list = "$file_i: $ln₍k₆₎ | $ln₍d₁₎ | $ln₍s₂₎\n"
end

open("Data/fakeData/params.txt", "a") do f 
    write(f, param_list)
end

println("Successfully created file at $(filename)")
