using Turing 
using Distributions
using Base.Threads
using ForwardDiff: Dual 
using ForwardDiff

include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/Differential/ode_restricted.jl")

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
        data, problem, selected_days, s, σ_likelihood, 
        log_norm::String,
        distro::ContinuousDistribution ; 
        num_experiments = 1)

    # Process the transform 
    @match log_norm begin 
        "identity" => global transform = (x -> x)
        "logarithm" => global transform = (x -> log(x))
    end
    # Process data matrix
    data = transform.(data)
    
    ## Regular priors
    ln_k6 ~ truncated(distro; lower=-100, upper=0)

    ## Convert ForwardDiff to Float64 (bad type interface)
    p = [ln_k6] .|> exp
    float_p = Vector{Float64}(undef, length(p))
    for (i, param) in enumerate(p) 
        float_p[i] = (typeof(param) <: Dual) ? param.value : param
    end

    ## Solve DDE model  
    repacked_p = repack_params(updateParams1(float_p...) )
    predictions = solve(problem; p=repacked_p[1], saveat=0.1)
    pred_vol = predictions[4,:] + predictions[5,:]
    sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]

    ## Likelihoods
    for exp in 1:num_experiments
        for i in eachindex(sliced_pred)
            data[exp, i] ~ (Normal(transform(sliced_pred[i]),  σ_likelihood))
        end
    end 
end

@model function fit_individual_restricted3(
    data, problem, selected_days, s, σ_likelihood, 
    log_norm::String,
    # distro::Vector{ContinuousDistribution}; 
    distro::ContinuousDistribution ;
    num_experiments = 1)

    # Process the transform 
    @match log_norm begin 
        "identity" => global forward = (x -> x)
        "logarithm" => global forward = (x -> log(x))
    end
    @match log_norm begin 
        "identity" => global reciprocal = (x -> x)
        "logarithm" => global reciprocal = (x -> exp(x))
    end
        
    ## Regular priors
    ln_k6 ~ truncated(distro; lower=-100, upper=0) # Negative half-Cauchy
    ln_d1 ~ truncated(distro; lower=0,    upper=7) # Positive half-Cauchy
    ln_s2 ~ truncated(distro; lower=-100, upper=0)

    ## Convert ForwardDiff to Float64 (bad type interface)
    p = [ln_k6, ln_d1, ln_s2] .|> exp
    float_p = Vector{Float64}(undef, length(p))
    for (i, param) in enumerate(p) 
        float_p[i] = (typeof(param) <: Dual) ? param.value : param
    end

    ## Solve DDE model  
    repacked_p = repack_params(updateParams3(float_p...) )
    predictions = solve(problem; p=repacked_p[1], saveat=0.1)
    pred_vol = predictions[4,:] + predictions[5,:]
    sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]

    ## Likelihoods
    for exp in 1:num_experiments
        for i in eachindex(sliced_pred)
            data[exp, i] ~ reciprocal(Normal(forward(sliced_pred[i]),  σ_likelihood))
        end
    end 
end
