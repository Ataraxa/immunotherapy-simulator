"""
Simple Bayesian model to fit the data to the simulated immune response.

It is NOT the hierarchical model, and is only used to fit one single 
time series, eg. for validation.

Assumptions that are being tested or will be in the future: (✓□✗)
    □ log-transformation of the search space
    ✓ truncation of hyperparameters vs regular parameters: not needed
"""

using Turing
using LinearAlgebra

# Custom libraries
include("ode_model.jl")
include("binormal.jl")

@model function fit_immune_resp(data, problem)
    if data === missing 
        data = Vector{Float64}(undef, 271)
    end

    # Prior distributions
    # k6
    mu1_k6 ~ Normal(0.2, 0.2)
    mu2_k6 ~ Normal(0.8, 0.2)
    s1_k6 ~ truncated(Normal(0, 0.2); lower=0)
    s2_k6 ~ truncated(Normal(0, 0.2); lower=0)
    a_k6 ~ Beta(1, 1)
    k6 ~ BiNormal{Float64}(mu1_k6, mu1_k6, s1_k6, s2_k6, a_k6)

    # d1 
    mu1_d1 ~ Normal(8, 2)
    mu2_d1 ~ Normal(14, 2)
    s1_d1 ~ truncated(Normal(0, 2); lower=0)
    s2_d1 ~ truncated(Normal(0, 2); lower=0)
    a_d1 ~ Beta(1, 1)
    d1 ~ BiNormal{Float64}(mu1_d1, mu1_d1, s1_d1, s2_d1, a_d1)
    
    # s2
    mu1_s2 ~ Normal(0.0, 0.2)
    mu2_s2 ~ Normal(0.6, 0.2)
    s1_s2 ~ truncated(Normal(0, 0.2); lower=0)
    s2_s2 ~ truncated(Normal(0, 0.2); lower=0)
    a_s2 ~ Beta(1, 1)
    s2 ~ BiNormal{Float64}(mu1_s2, mu1_s2, s1_s2, s2_s2, a_s2)
    
    # Simulate immune response
    p=[k6, d1, s2]
    predictions = solve(problem, MethodOfSteps(Tsit5()); p=p, saveat=0.1)
    tumour_vol_prediction = predictions[4,:] + predictions[5,:]
    # display(plot(tumour_vol_prediction))
    
    # Experimental error (σ_err)
    σ_err = 0.1 

    # Observations.
    for i in eachindex(tumour_vol_prediction)
        data[i] ~ Normal(tumour_vol_prediction[i], σ_err^2)
    end

end