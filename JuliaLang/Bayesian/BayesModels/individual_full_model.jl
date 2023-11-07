using Turing 
using Distributions 
using ForwardDiff: Dual 
using ForwardDiff

include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")

"""
Statistical model to estimate all the individual parameters of a mouse, given an
experimental time series (θ from R21, ie not including the initial 
conditions u0).

Inputs:

    - problem: a tumour DDEProblem that requires a full parameter vector (w/o 
    u0) to be solved
    - data: the time series to be fitted
    - selected_days days: vector of Int64, days at which data was collected
    - timestep: to solve the DDEProblem and synchronise with selected_days
    - σ_likelihood: standard deviation of the likelihood distribution
"""
@model function fit_individual(data, problem, selected_days, timestep, 
        σ_likelihood)

    # Prior distributions
    c = christian

    ln_td      ~ truncated(Normal(c.t_d,     1); lower=-100, upper=log(2000))
    ln_t_delay ~ truncated(Normal(c.t_delay, 1); lower=-100, upper=log(2000))
    ln_t_last  ~ truncated(Normal(c.t_last,  1); lower=-100, upper=log(2000))
   
    ln_t_delay12 = 0 # Not useful for current treatment
    ln_t_last12 = 0

    ln_k1 ~ truncated(Normal(c.k1,1); lower=-100, upper=log(2000))
    ln_k2 ~ truncated(Normal(c.k2,1); lower=-100, upper=log(2000))
    ln_k3 ~ truncated(Normal(c.k3,1); lower=-100, upper=log(2000))
    ln_k4 ~ truncated(Normal(c.k4,1); lower=-100, upper=log(2000))
    ln_k5 ~ truncated(Normal(c.k5,1); lower=-100, upper=log(2000))
    ln_k6 ~ truncated(Normal(c.k6,1); lower=-100, upper=log(2000))
    
    ln_d1 ~ truncated(Normal(c.d1,1); lower=-100, upper=log(2000))
    ln_d2 ~ truncated(Normal(c.d2,1); lower=-100, upper=log(2000))
    ln_d3 ~ truncated(Normal(c.d3,1); lower=-100, upper=log(2000))
    ln_d4 ~ truncated(Normal(c.d4,1); lower=-100, upper=log(2000))
    ln_d5 ~ truncated(Normal(c.d5,1); lower=-100, upper=log(2000))
    ln_d6 ~ truncated(Normal(c.d6,1); lower=-100, upper=log(2000))
    ln_d7 ~ truncated(Normal(c.d7,1); lower=-100, upper=log(2000))
    ln_d8 ~ truncated(Normal(c.d8,1); lower=-100, upper=log(2000))

    ln_s1 ~ truncated(Normal(c.s1,1); lower=-100, upper=log(2000))
    ln_s2 ~ truncated(Normal(c.s2,1); lower=-100, upper=log(2000))

    # DDEProblem behaves badly with ForwardDiff, so need to check types 
    p = [ln_td, ln_t_delay, ln_t_last, ln_t_delay12, ln_t_last12,
        ln_k1, ln_k2, ln_k3, ln_k4, ln_k5, ln_k6, 
        ln_d1, ln_d2, ln_d3, ln_d4, ln_d5, ln_d6, ln_d7, ln_d8, 
        ln_s1, ln_s2]

    p = exp.(p) # Element-wise exponentiation of the array
    float_p = Vector{Float64}(undef, length(p))
    for (i, param) in enumerate(p)
        if typeof(param) <: Dual
            float_p[i] = param.value
        end
    end

    # Solve the DDE model
    predictions = solve(problem; p=float_p, saveat=0.1)
    pred_tumour = predictions[4,:] + predictions[5,:]
    println("$(size(predictions)) | $(size(pred_tumour))")
    if size(pred_tumour) == (1,)
        print(float_p)
    end
    # println(pred_tumour)
    sliced_pred = pred_tumour[selected_days*trunc(Int, 1/timestep) .+ 1]
    # sliced_pred = pred_tumour[[1, 71, 81, 91, 111, 141, 171, 201]]

    # Likelihoods
    for i in eachindex(sliced_pred)
        data[1, i] ~ Normal(log(sliced_pred[i]), σ_likelihood)
    end    
end