include("../Model/Differential/ode_core.jl")

prob = problem_factory()
solve(prob)

println("Test was passed successfully!")