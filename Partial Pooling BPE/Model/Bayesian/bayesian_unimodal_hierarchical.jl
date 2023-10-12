using Turing 
using Distributions
using LinearAlgebra
using Base.Threads
using ForwardDiff: Dual 
using ForwardDiff
include("../Differential/ode_model.jl")

"""
Statistical hierachical model: population prior is unimodal (with 2 hyperparameters)

Important: the `num_experiments` parameter can be used to slice down the data matrix,
so that only the first _n_ rows are used.
"""
@model function fit_unimodal_hierarchical(data, problem, num_experiments, s, selected_days)

    # println("Starting evaluation of the model: ", threadid())

    # Initialise the parameter arrays
    k6 = Vector{Float64}(undef, num_experiments)
    d1 = Vector{Float64}(undef, num_experiments)
    s2 = Vector{Float64}(undef, num_experiments)

    # Hyperprior distributions
    ln_µ_k6 ~ Normal(-0.70, 1)
    ln_σ_k6 ~ truncated(Normal(0, 1); lower=0)

    ln_µ_d1 ~ Normal(2.40, 1)
    ln_σ_d1 ~ truncated(Normal(0, 1); lower=0)

    ln_µ_s2 ~ Normal(-1.2, 1)
    ln_σ_s2 ~ truncated(Normal(0, 1); lower=0)
    
    # Regular priors
    for expr in 1:num_experiments
        k6[expr] ~ truncated(Normal(ln_µ_k6, ln_σ_k6); lower=-100, upper=log(10))
        d1[expr] ~ truncated(Normal(ln_µ_d1, ln_σ_d1); lower=-100, upper=log(75))
        s2[expr] ~ truncated(Normal(ln_µ_s2, ln_σ_s2); lower=-100, upper=log(10))    
    end

    # Experimental error (σ_err)
    σ_err = 10

    # Likelihood 
    # println("Gonna go into likelihoods")
    # expr = experiment n°X
    for expr in 1:num_experiments
        p = [exp(k6[expr]), exp(d1[expr]), exp(s2[expr])]
        valued_p = Vector{Float64}(undef, 3)
        must_switch = false
        
        # Must convert FowardDiff to Float numbers
        # DDESolver behaves badly with ForwardDiff type
        for (i, param) in enumerate(p) 
            if typeof(param) == Dual{ForwardDiff.Tag{Turing.TuringTag, Float64}, Float64, 3}
                valued_p[i] = param.value 
                println("Changed!")
                must_switch = true
            end 
        end
    
        predictions = solve(problem, MethodOfSteps(Tsit5()); p=(must_switch) ? valued_p : p, saveat=0.1)
        pred_vol = predictions[4,:] + predictions[5,:]
        sliced_pred = pred_vol[selected_days*trunc(Int, 1/s) .+ 1]

        for i in eachindex(sliced_pred)
            # println("One timestep...")
            data[expr, i] ~ Normal(sliced_pred[i], σ_err)
        end
    end

    # println("Finished evaluation")
end
