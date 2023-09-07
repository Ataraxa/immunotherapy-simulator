"""
This files simply contains a dummy problem, created for a Bayesian inference tutorial 
(see: https://turing.ml/v0.22/tutorials/10-bayesian-differential-equations/)

This is useful to check if all the modules work correctly, and to test some 
hypothesis on a small-scale model.

Running the full inference takes approx. 1 min on a normal laptop.
"""

using Turing
using DifferentialEquations
using StatsPlots: plot
using LinearAlgebra

# Set a seed for reproducibility
using Random
Random.seed!(14);

# DDE Model
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
    σ ~ InverseGamma(2, 3)
    α ~ truncated(Normal(1.5, 0.5); lower=0.5, upper=2.5)
    β ~ truncated(Normal(1.2, 0.5); lower=0, upper=2)
    γ ~ truncated(Normal(3.0, 0.5); lower=1, upper=4)
    δ ~ truncated(Normal(1.0, 0.5); lower=0, upper=2)

    # Simulate Lotka-Volterra model.
    p = [α, β, γ, δ]
    predicted = solve(prob, MethodOfSteps(Tsit5()); p=p, saveat=0.1)

    # Observations.
    for i in eachindex(predicted)
        data[:, i] ~ MvNormal(predicted[i], σ^2 * I)
    end
end

# Bayesian inference
model_dde = fitlv_dde(ddedata, prob_dde)

# Sample 3 independent chains
chain_dde = sample(model_dde, NUTS(0.65), MCMCSerial(), 300, 3; progress=false)

plot(chain_dde)

