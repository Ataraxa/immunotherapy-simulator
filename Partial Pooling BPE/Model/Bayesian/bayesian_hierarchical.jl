using Turing 

include("ode_model.jl")
include("binormal.jl")

@model function fit_hierarchical(data, problem, num_experiments)
    # Initialise the parameter arrays
    k6 = Vector{Float64}(undef, num_experiments)
    d1 = Vector{Float64}(undef, num_experiments)
    s2 = Vector{Float64}(undef, num_experiments)

    # Hyperprior distributions
    µ1_k6 ~ Normal()
    µ2_k6 ~ Normal()
    σ1_k6 ~ truncated(Normal())
    σ2_k6 ~ truncated(Normal())
    α_k6 ~ Beta(1, 1)

    µ1_d1 ~ Normal()
    µ2_d1 ~ Normal()
    σ1_d1 ~ truncated(Normal())
    σ2_d1 ~ truncated(Normal())
    α_d1 ~ Beta(1, 1)

    µ1_s2 ~ Normal()
    µ2_s2 ~ Normal()
    σ1_s2 ~ truncated(Normal())
    σ2_s2 ~ truncated(Normal())
    α_s2 ~ Beta(1, 1)
    
    # Regular priors
    for exp in data
        k6[exp] ~ BiNormal(µ1_k6, µ2_k6, σ1_k6, σ2_k6, α_k6)
        d1[exp] ~ BiNormal(µ1_d1, µ2_d1, σ1_d1, σ2_d1, α_d1)
        s2[exp] ~ BiNormal(µ1_s2, µ2_s2, σ1_s2, σ2_s2, α_s2)
    end

    # Experimental error (σ_err)
    σ_err = 0.1 

    # Likelihood 
    for exp in data 
        p = [k6[exp], d1[exp], s2[exp]]
        predictions = solve(problem, MethodOfSteps(Tsit5()); p=p, saveat=0.1)
        pred_vol = predictions[4,:] + predictions[5,:]

        for i in eachindex(pred_vol)
            data[exp, i] ~ Normal(pred_vol[i], σ_err^2)
        end
    end
end
