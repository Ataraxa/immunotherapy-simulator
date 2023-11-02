using Turing 
using Distributions
using Base.Threads
using ForwardDiff: Dual 
using ForwardDiff

include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/Differential/ode_restricted.jl")

"""
Statistical hierachical model: population prior is unimodal (with 2 hyperparameters)

Important: the `num_experiments` parameter can be used to slice down the data matrix,
so that only the first _n_ rows are used.
"""
@model function fit_individual_restricted(data, problem, s, selected_days, 
        std_k6, std_d1, std_s2,
        max1, max2,
        exp_err)
    
    # println("__________________________________________________")

    ## Regular priors
    ln_k6 ~ truncated(Normal(-0.7, std_k6); lower=-100, upper=log(max1))
    ln_d1 ~ truncated(Normal( 2.3, std_d1); lower=-100, upper=log(max2))
    ln_s2 ~ truncated(Normal(-1.3, std_s2); lower=-100, upper=log(max1))    

    ## Experimental error (σ_err)
    σ_err = exp_err

    ## Convert ForwardDiff to Float64 (bad type interface)
    p = [ln_k6, ln_d1, ln_s2] .|> exp
    float_p = Vector{Float64}(undef, length(p))
    for (i, param) in enumerate(p) 
        float_p[i] = (typeof(param) <: Dual) ? param.value : param
    end

    ## Solve DDE model  
    update = updateParams(float_p...)
    problem = restricted_simulation(update)  
    predictions = solve(problem; saveat=0.1)
    pred_vol = predictions[4,:] + predictions[5,:]
    sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]

    ## Likelihoods
    for i in eachindex(sliced_pred)
        data[1, i] ~ Normal(log(sliced_pred[i]), σ_err) # TODO: should it be log??
    end
end
