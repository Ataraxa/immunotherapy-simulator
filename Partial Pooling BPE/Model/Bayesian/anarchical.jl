using Turing 
using Distributions
using Base.Threads
include("../Differential/ode_model.jl")

"""
Statistical hierachical model: population prior is unimodal (with 2 hyperparameters)

Important: the `num_experiments` parameter can be used to slice down the data matrix,
so that only the first _n_ rows are used.
"""
@model function fit_anarchical(data, problem, s, selected_days, 
        std_k6, std_d1, std_s2,
        max1, max2,
        exp_err)
    
    printn("__________________________________________________")

    # Regular priors
    k6 ~ truncated(Normal(-0.7, std_k6); lower=-100, upper=log(max1))
    d1 ~ truncated(Normal( 2.3, std_d1); lower=-100, upper=log(max2))
    s2 ~ truncated(Normal(-1.3, std_s2); lower=-100, upper=log(max1))    

    # Experimental error (σ_err)
    σ_err = exp_err

    # Simulate model
    p = [exp(k6), exp(d1), exp(s2)]
    println(p)

    # Experimental block 
    for param in p 
        if typeof(param) == Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}, Float64, 3}
            param = param.value 
        end 
    end

    predictions = solve(problem; p=p, saveat=0.1)
    pred_vol = predictions[4,:] + predictions[5,:]
    sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]
    println(p)
    println(sliced_pred)
    for i in eachindex(sliced_pred)
        data[i] ~ Normal(log(sliced_pred[i]), σ_err)
    end
end
