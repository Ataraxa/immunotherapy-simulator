using Distributions

include("../../Model/Differential/ode_restricted.jl")


struct updateParams 
    k6
    d1
    s2
end

function do_simulations(n_iters=1_000)
    for i in 1:n_iters
        k6 = rand(Uniform(0.0, 1.0 ))
        d1 = rand(Uniform(7.0, 15.0))
        s2 = rand(Uniform(0.0, 1.0 ))

        update = updateParams(k6, d1, s2)
        prob = restricted_simulation(update)
        sol = solve(prob; saveat=0.1)
    end
end

do_simulations()