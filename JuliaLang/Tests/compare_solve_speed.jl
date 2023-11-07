using DifferentialEquations

include("../Model/Differential/ode_core.jl")
include("../Model/Differential/ode_params.jl")
include("../Model/Differential/ode_restricted.jl")

prob_dde = create_problem()
@time begin
    println("Unrestricted parameter space")
    sol = solve(prob_dde; saveat=0.1)
end

@time begin
    println("Restricted parameter space")
    p, u0 = repack_params(updateParams(0.1, 12, 0.3))
    sol = solve(prob_dde; p=p)
end

print("Test done")