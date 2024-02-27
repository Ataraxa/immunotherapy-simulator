#=
Takes a reference matrix with no input noise.
Outputs the same matrix but with added noise, and saves it to a file.
If needed, can add several repetitions, resulting in a tensor at the output. 
=#

using JLD

function add_noise(ref, rep::Int, noise_level::Float64)
    output = Array{Float64}(undef, size(ref)..., rep)

    for i in 1:rep
        noise_matrix = rand(Normal(0, 1)*noise_level, size(ref)...)
        output[:,:,i] = ref + noise_matrix 
    end

    return output
end

path = "Data/fake_data"

# Load ref data 
input = 3
ref = load("$path/trajectories-$input.jld", "M")[:,:,1]

# Add noise 
output = add_noise(ref, 10, 1.)

# Save file 
file_i = 0
filename = "trajectories-$file_i.jld"
while isfile("$path/$filename")
    global file_i+=1
    global filename = "trajectories-$file_i.jld"
end
save("$path/$filename", "M", output)

summary = "$(file_i) -> Copy of $input (with noise=1) \n"
open("$path/log.txt", "a") do f 
    write(f, summary)
end

open("$path/params.txt", "a") do f 
    write(f, "$(file_i) -> copy of $input")
end

println("Successfully created file at $(filename)")
