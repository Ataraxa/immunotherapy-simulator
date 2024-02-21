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
    
    ## Regular priors
    ln_k6 ~ distro[1] # Negative half-Cauchy
    ln_d1 ~ distro[2] # Positive half-Cauchy
    ln_s2 ~ distro[3]

    ## Convert ForwardDiff to Float64 (bad type interface)
    p = [ln_k6, ln_d1, ln_s2] .|> exp
    # float_p = p

    float_p = Vector{Float64}(undef, length(p))
    for (i, param) in enumerate(p) 
        if typeof(param) <: Dual
            float_p[i] = param.value + sum([i for i in param.partials])
        else
            float_p[i] = param
        end
    end

    ## Solve DDE model  
    params[var_params_index] .= float_p
    pred = solve(problem; p=params, saveat=s)

    v = sum(pred[4:end,:], dims=1)
    combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
    if size(combined_pred, 2) != 271
        println.(p)
        println("______________________")
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
