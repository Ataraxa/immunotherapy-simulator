using Turing 

include("../Differential/ode_model.jl")
include("../binormal.jl")

@model function fit_hierarchical(data, problem, num_experiments, s, selected_days)
    # Initialise the parameter arrays
    k6 = Vector{Float64}(undef, num_experiments)
    d1 = Vector{Float64}(undef, num_experiments)
    s2 = Vector{Float64}(undef, num_experiments)

    # Hyperprior distributions
    µ1_k6 ~ Normal(log(0.2), 1)
    µ2_k6 ~ Normal(log(0.5), 1)
    σ1_k6 ~ truncated(Normal(0, 1); lower=0)
    σ2_k6 ~ truncated(Normal(0, 1); lower=0)
    α_k6 ~ Beta(1, 1)

    µ1_d1 ~ Normal(log( 7), 1)
    µ2_d1 ~ Normal(log(15), 1)
    σ1_d1 ~ truncated(Normal(0, 1); lower=0)
    σ2_d1 ~ truncated(Normal(0, 1); lower=0)
    α_d1 ~ Beta(1, 1)

    µ1_s2 ~ Normal(log(0.2), 1)
    µ2_s2 ~ Normal(log(0.5), 1)
    σ1_s2 ~ truncated(Normal(0, 1); lower=0)
    σ2_s2 ~ truncated(Normal(0, 1); lower=0)
    α_s2 ~ Beta(1, 1)
    
    # Regular priors
    for expr in data
        k6[expr] ~ BiNormal(µ1_k6, µ2_k6, σ1_k6, σ2_k6, α_k6)
        d1[expr] ~ BiNormal(µ1_d1, µ2_d1, σ1_d1, σ2_d1, α_d1)
        s2[expr] ~ BiNormal(µ1_s2, µ2_s2, σ1_s2, σ2_s2, α_s2)
    end

    # Experimental error (σ_err)
    σ_err = 4

    # Likelihood 
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
            data[expr, i] ~ Normal(sliced_pred[i], 4)
        end
    end
end
