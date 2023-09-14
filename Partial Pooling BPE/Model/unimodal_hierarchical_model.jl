using Turing 
using Distributions
include("ode_model.jl")
include("binormal.jl")

"Dummy statistical hierachical model"
@model function fit_dummy_hierarchical(data, problem, num_experiments, s, selected_days)
    # Initialise the parameter arrays
    k6 = Vector{Float64}(undef, num_experiments)
    d1 = Vector{Float64}(undef, num_experiments)
    s2 = Vector{Float64}(undef, num_experiments)

    # Hyperprior distributions
    µ_k6 ~ Normal(0.5, 0.2)
    σ_k6 ~ truncated(Normal(0, 0.3); lower=0)

    µ_d1 ~ Normal(11, 2)
    σ_d1 ~ truncated(Normal(0, 3); lower=0)

    µ_s2 ~ Normal(0.3, 0.2)
    σ_s2 ~ truncated(Normal(0, 0.2); lower=0)
    
    # Regular priors
    for exp in data
        k6[exp] ~ truncated(Normal(µ_k6, σ_k6); lower=0)
        d1[exp] ~ truncated(Normal(µ_d1, σ_d1); lower=0)
        s2[exp] ~ truncated(Normal(µ_s2, σ_s2); lower=0)    
    end

    # Experimental error (σ_err)
    σ_err = 0.1 

    # Likelihood 
    for exp in data 
        p = [k6[exp], d1[exp], s2[exp]]
        predictions = solve(problem, MethodOfSteps(Tsit5()); p=p, saveat=0.1)
        pred_vol = predictions[4,:] + predictions[5,:]
        sliced_pred = pred_vol[selected_days*1/s .+ 1]

        for i in eachindex(sliced_pred)
            data[exp, i] ~ Normal(sliced_pred[i], σ_err^2)
        end
    end
end
