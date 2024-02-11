using Turing 
using Distributions
using Base.Threads
using ForwardDiff: Dual 
using ForwardDiff

include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")

"""
Statistical model to estimate only the sensitive parameters of a mouse, given an
experimental time series (θ from ~R3, ie not including the initial 
conditions u0).

Inputs:

    - problem: a tumour DDEProblem that requires a full parameter vector (w/o 
    u0) to be solved
    - data: the time series to be fitted
    - selected_days days: vector of Int64, days at which data was collected
    - timestep: to solve the DDEProblem and synchronise with selected_days
    - σ_likelihood: standard deviation of the likelihood distribution
"""
@model function fit_individual_restricted1(
        data, problem, selected_days, s, σ_err, 
        log_norm::SubString{String},
        distro::Vector{T} where T <: ContinuousDistribution ; 
        num_experiments = 1)

    # Process the transform 
    @match log_norm begin 
        "norm_noise" => global transform = (x -> x)
        "logn_noise" => global transform = (x -> log(x))
    end
    # Process data matrix
    data = transform.(data)
    
    ## Regular priors
    ln_k6 ~ truncated(distro[1]; lower=-100, upper=0)

    ## Convert ForwardDiff to Float64 (bad type interface)
    p = [ln_k6] .|> exp
    float_p = Vector{Float64}(undef, length(p))
    for (i, param) in enumerate(p) 
        float_p[i] = (typeof(param) <: Dual) ? param.value : param
    end

    ## Solve DDE model  
    params = copy(christian_true_params)
    params
    predictions = solve(problem; p=repacked_p[1], saveat=0.1)
    pred_vol = predictions[4,:] + predictions[5,:]
    sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]

    ## Likelihoods
    for exp in 1:num_experiments
        for i in eachindex(sliced_pred)
            data[exp, i] ~ (Normal(transform(sliced_pred[i]),  σ_err))
        end
    end 
end

@model function fit_individual_restricted3(
    data, problem, selected_days, s, σ_err, 
    log_norm::SubString{String},
    distro::Vector{T} where T <: ContinuousDistribution,
    var_params_index; 
    num_experiments = 1)

    # Process the transform 
    @match log_norm begin 
        "none" => global transform = (x -> x)
        "loga" => global transform = (x -> log(x)) # logan paul?
    end
    # Process data matrix
    data = transform.(data)
    params = copy(christian_true_params) 
    
    # While-loop to avoi stiff paramert combinations
    pred = Array{Float64}(undef, 5, 10)
    while size(pred)[2] != 271
        ## Regular priors
        ln_k6 ~ truncated(distro[1]; lower=-100, upper=7) # Negative half-Cauchy
        ln_d1 ~ truncated(distro[2]; lower=-5,   upper=7) # Positive half-Cauchy
        ln_s2 ~ truncated(distro[3]; lower=-100, upper=7)

        ## Convert ForwardDiff to Float64 (bad type interface)
        p = [ln_k6, ln_d1, ln_s2] .|> exp
        float_p = Vector{Float64}(undef, length(p))
        for (i, param) in enumerate(p) 
            float_p[i] = (typeof(param) <: Dual) ? param.value : param
        end
        ## Solve DDE model  
        
        params[var_params_index] .= float_p
    
        pred = solve(problem; p=params, saveat=s)
        if size(pred)[2] != 271
            println.(float_p)
            println("___________________________")
        end
    end
    v = sum(pred[4:end,:], dims=1)
    combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
    if size(combined_pred, 2) != 271
        println(float_p)
    end
    sliced_pred = combined_pred[:,selected_days*trunc(Int, 1/s) .+ 1]
    
    ## Likelihood
    for exp in 1:num_experiments
        for qty in 1:4 # range g > c > p > v
            for i in eachindex(sliced_pred[1,:])
                data[qty,i,exp] ~ Normal(transform(sliced_pred[qty,i]),  σ_err)
            end
        end
    end 
end
