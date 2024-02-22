# Define Lotka-Volterra model.
function lotka_volterra(du, u, p, t)
    # Model parameters.
    α, β, γ, δ = p
    # Current state.
    x, y = u

    # Evaluate differential equations.
    du[1] = (α - β * y) * x # prey
    du[2] = (δ * x - γ) * y # predator

    return nothing
end

# Define initial-value problem.
u0 = [1.0, 1.0]
p = [1.5, 1.0, 3.0, 1.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lotka_volterra, u0, tspan, p)

# Define noise level
σ = 0.8
sol = solve(prob, Tsit5(); saveat=0.1)
odedata = Array(sol) + σ * randn(size(Array(sol)))

# Plot simulation and noisy observations.
# plot(sol; alpha=0.3)
# scatter!(sol.t, odedata'; color=[1 2], label="")

### Save predictions to file
file_i = 0
path = "Data/fake_data"
filename = "trajectories-$file_i.jld"
while isfile("$path/$filename")
    global file_i+=1
    global filename = "trajectories-$file_i.jld"
end
# if do_overwrite_prev
#     filename = "trajectories-$(file_i-1).jld"
#     global file_i -= 1
# end
save("$path/$filename", "M", odedata)

summary = "$(file_i) ->  N/A | $σ | +(.) | Lotka-Volterra \n"
open("$path/log.txt", "a") do f 
    write(f, summary)
end

open("$path/params.txt", "a") do f 
    write(f, "$(file_i) -> "*join(string.([1.5, 1.0, 3.0, 1.0]), " | ")*"\n")
end

println("Successfully created file at $(filename)")