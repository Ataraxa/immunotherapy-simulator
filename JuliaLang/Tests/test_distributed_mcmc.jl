using Turing
using DifferentialEquations
using StatsPlots: plot
using LinearAlgebra
using Random
Random.seed!(14);

# DDE Model
function parallel_wrapper(invgamma, normal, trunca, ts5, solve, MethodOfSteps, I, MvNormal)

    function delay_lotka_volterra(du, u, h, p, t)
        # Model parameters.
        α, β, γ, δ = p

        # Current state.
        x, y = u
        # Evaluate differential equations
        du[1] = α * h(p, t - 1; idxs=1) - β * x * y
        du[2] = -γ * y + δ * x * y

        return nothing
    end

    # Define initial-value problem
    p = (1.5, 1.0, 3.0, 1.0)
    u0 = [1.0; 1.0]
    tspan = (0.0, 10.0)
    h(p, t; idxs::Int) = 1.0
    prob_dde = DDEProblem(delay_lotka_volterra, u0, h, tspan, p);

    # Solve the problem and generate a random approxiation of the solution
    sol_dde = solve(prob_dde; saveat=0.1)
    ddedata = Array(sol_dde) + 0.5 * randn(size(sol_dde))

    # Plot simulation and noisy observations
    # StatsPlots.plot(sol_dde)
    # scatter!(sol_dde.t, ddedata'; color=[1 2], label="")

    # Statistical model
    @model function fitlv_dde(data, prob)

        # Prior distributions.
        σ ~ invgamma(2, 3)
        α ~ trunca(normal(1.5, 0.5); lower=0.5, upper=2.5)
        β ~ trunca(normal(1.2, 0.5); lower=0, upper=2)
        γ ~ trunca(normal(3.0, 0.5); lower=1, upper=4)
        δ ~ trunca(normal(1.0, 0.5); lower=0, upper=2)

        # Simulate Lotka-Volterra model.
        p = [α, β, γ, δ]
        predicted = solve(prob, MethodOfSteps(ts5()); p=p, saveat=0.1)

        # Observations.
        for i in eachindex(predicted)
            data[:, i] ~ MvNormal(predicted[i], σ^2 * I)
        end
    end

    # Bayesian inference
    model_dde = fitlv_dde(ddedata, prob_dde)

    # Test sampling
    # display(@btime chain_dde_nuts = sample(model_dde, NUTS(0.65), 100; progress=false))
    chain_dde_nuts = sample(fitlv_dde(ddedata, prob_dde), NUTS(0.65), MCMCDistributed(),
        100, 1; progress=false)
end

parallel_wrapper(InverseGamma, Normal, truncated, Tsit5, solve, MethodOfSteps, I, MvNormal)