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
    base_params::Vector{Float64},
    log_norm::SubString{String},
    distro::Vector{T} where T <: ContinuousDistribution,
    var_params_index; 
    num_experiments = 1)

    # Process the transform [✓]
    @match log_norm begin 
        "none" => global transform = (x -> x)
        "loga" => global transform = (x -> log(x)) # logan paul?
    end
    data = transform.(data)

    ## Baseline parameters [✓]
    params = Vector{Any}(undef, length(base_params))
    params[:] = copy(base_params) 
    
    ## Regular priors [✓]
    ln_k6 ~ distro[1]
    ln_d1 ~ distro[2]
    ln_s2 ~ distro[3]
        
    ## Solve DDE model [no]
    p = [ln_k6, ln_d1, ln_s2] .|> exp
    params[var_params_index] .= p
    pred = solve(problem; p=params, saveat=s)
    # pred = solve(problem, AutoTsit5(RadauIIA3(); 
    #                         maxstiffstep=70, stifftol=1.4, # low stifftol -> everything is stiff
    #                         maxnonstiffstep=1, nonstifftol=0.7, 
    #                         stiffalgfirst=true); 
    #                 p=params, saveat=s)
    if size(pred, 2) != 271
        println.(p)
        println("______________________")
        global pred = solve(problem, RadauIIA5(); p=params, saveat=s)
    end

    ## Process tumour data (sum + slice at active days) [no]
    v = sum(pred[4:end,:], dims=1)
    combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
    sliced_pred = combined_pred[:,selected_days*trunc(Int, 1/s) .+ 1]
    
    ## Likelihood [✓]
    for exp in 1:num_experiments
        for qty in 1:4 # range g > c > p > v
            for i in eachindex(sliced_pred[1,:])
                data[qty,i,exp] ~ Normal(transform(sliced_pred[qty,i]),  σ_err)
            end
        end
    end 
end

@model function fit_lotkaVolterra(
    data, problem, selected_days, stepsize, σ_err, 
    true_params::Vector{Float64},
    log_norm::SubString{String}, # To transform the data
    distro::Vector{T} where T <: ContinuousDistribution, # Vector of priors
    var_params_index; 
    num_experiments = 1
    )

    # Process the transform 
    @match log_norm begin 
        "none" => global transform = (x -> x)
        "+(.)" => global transform = (x -> x)
        "loga" => global transform = (x -> log(x)) # logan paul?
    end
    data = transform.(data)

    ## Baseline parameters
    params = Vector{Any}(undef, length(true_params))
    params[:] = copy(true_params) 
    
    ## Regular priors
    α ~ distro[1]
    β ~ distro[2]
    γ ~ distro[3]
    δ ~ distro[4]

    # ## Convert ForwardDiff to Float64 (bad type interface)
    # p = [ln_k6, ln_d1, ln_s2] .|> exp

    ## Solve DDE model  
    params[var_params_index] .= [α, β, γ, δ]
    pred = solve(problem; p=params, saveat=stepsize)
    sliced_pred = pred[:,[1, 20, 30, 40, 50, 60, 70, 80]]
    ## Process tumour data (sum + slice at active days)
    # v = sum(pred[4:end,:], dims=1)
    # combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
    # if size(combined_pred, 2) != 271
    #     println.(p)
    #     println("______________________")
    # end
    # sliced_pred = combined_pred[:,selected_days*trunc(Int, 1/s) .+ 1]
    
    ## Likelihood
    for exp in 1:num_experiments
        for qty in 1:2 # range g > c > p > v TODO:change back to 4
            for i in eachindex(sliced_pred[1,:])
                data[qty,i,exp] ~ Normal(transform(sliced_pred[qty,i]),  σ_err)
            end
        end
    end 
end