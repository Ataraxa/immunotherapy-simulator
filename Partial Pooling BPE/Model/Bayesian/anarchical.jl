using Turing 
using Distributions
using Base.Threads
include("../Differential/ode_model.jl")

"""
Statistical hierachical model: population prior is unimodal (with 2 hyperparameters)

Important: the `num_experiments` parameter can be used to slice down the data matrix,
so that only the first _n_ rows are used.
"""
@model function fit_anarchical(data, problem, s, selected_days, upper1, upper2)
    # Regular priors
    d1 ~ truncated(Normal(2.3, 0.6);  lower=-100, upper=log(5))
    s2 ~ truncated(Normal(-1.3, 0.6); lower=-100, upper=log(75))    
    k6 ~ truncated(Normal(-0.7, 0.6); lower=-100, upper=log(5))

    # Experimental error (σ_err)
    σ_err = 10 

    # Simulate model
    p = [exp(k6), exp(d1), exp(s2)]

    # sleep(3)
    # println(p)
    # println("$k6, $d1, $s2")
    # println("__________________")
    
    predictions = solve(problem; p=p, saveat=0.1)
    pred_vol = predictions[4,:] + predictions[5,:]
    sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]

    for i in eachindex(sliced_pred)
        data[i] ~ Normal(sliced_pred[i], σ_err^2)
    end
end
