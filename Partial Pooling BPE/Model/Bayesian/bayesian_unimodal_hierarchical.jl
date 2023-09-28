using Turing 
using Distributions
using LinearAlgebra
using Base.Threads
include("../Differential/ode_model.jl")

"""
Statistical hierachical model: population prior is unimodal (with 2 hyperparameters)

Important: the `num_experiments` parameter can be used to slice down the data matrix,
so that only the first _n_ rows are used.
"""
@model function fit_unimodal_hierarchical(data, problem, num_experiments, s, selected_days)
    println("Starting evaluation of the model: ", threadid())

    # Initialise the parameter arrays
    k6 = Vector{Float64}(undef, num_experiments)
    d1 = Vector{Float64}(undef, num_experiments)
    s2 = Vector{Float64}(undef, num_experiments)

    # Hyperprior distributions
    µ_k6 ~ Normal(-0.70, 1)
    σ_k6 ~ truncated(Normal(0, 1); lower=0)

    µ_d1 ~ Normal(2.40, 1)
    σ_d1 ~ truncated(Normal(0, 1); lower=0)

    µ_s2 ~ Normal(-1.2, 1)
    σ_s2 ~ truncated(Normal(0, 1); lower=0)
    
    # Regular priors
    for expr in 1:num_experiments
        k6[expr] ~ truncated(Normal(µ_k6, σ_k6); lower=-100, upper=log(10))
        d1[expr] ~ truncated(Normal(µ_d1, σ_d1); lower=-100, upper=log(75))
        s2[expr] ~ truncated(Normal(µ_s2, σ_s2); lower=-100, upper=log(10))    
    end

    # Experimental error (σ_err)
    σ_err = 10

    # Likelihood 
    # println("Gonna go into likelihoods")
    for expr in 1:num_experiments
        p = [exp(k6[expr]), exp(d1[expr]), exp(s2[expr])]
        predictions = solve(problem, MethodOfSteps(Tsit5()); p=p, saveat=0.1)
        pred_vol = predictions[4,:] + predictions[5,:]
        sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]

        for i in eachindex(sliced_pred)
            # println("One timestep...")
            data[expr, i] ~ Normal(sliced_pred[i], σ_err)
        end
    end

    # println("Finished evaluation")
end
