using Turing
using LinearAlgebra

# Custom libraries
include("ode_model.jl")
include("binormal.jl")

@model function fit_immune_resp(data, problem)
    # Prior distributions
    # k6
    mu1_k6 ~ Normal(0.2, 0.1)
    mu1_k6 ~ Normal(0.8, 0.1)
    s1_k6 ~ truncated(Normal(0, 0.1); lower=0)
    s2_k6 ~ truncated(Normal(0, 0.1); lower=0)
    a_k6 ~ Beta(1, 1)
    k6 ~ BiNormal{Float64}(mu1_k6, mu1_k6, s1_k6, s2_k6, a_k6)

    # d1 
    mu1_d1 ~ Normal(0.2, 0.1)
    mu2_d1 ~ Normal(0.8, 0.1)
    s1_d1 ~ truncated(Normal(0, 0.1); lower=0)
    s2_d1 ~ truncated(Normal(0, 0.1); lower=0)
    d1 ~ BiNormal{Float64}(mu1_d1, mu1_d1, s1_d1, s2_d1, a_d1)
    
    # s2
    mu1_s2 ~ Normal(0.2, 0.1)
    mu2_s2 ~ Normal(0.8, 0.1)
    s1_s2 ~ truncated(Normal(0, 0.1); lower=0)
    s2_s2 ~ truncated(Normal(0, 0.1); lower=0)
    s2 ~ BiNormal{Float64}(mu1_s2, mu1_s2, s1_s2, s2_s2, a_s2)

    # Experimental error
    s_err ~ truncated(Normal(0, 10); lower=0)
    
    # Simulate immune response
    _, p = get_default_values()
    p[2][6] = k6
    p[3][1] = d1 
    p[4][2] = s2
    predictions = solve(problem, MethodOfSteps(Tsit5()); p=p, saveat=0.1)
    tumour_vol_prediction = predictions[4] + predictions[5]

    # Observations.
    for i in 1:length(predicted)
        data[:, i] ~ Normal(tumour_vol_prediction[i], σ^2 * I)
    end
end