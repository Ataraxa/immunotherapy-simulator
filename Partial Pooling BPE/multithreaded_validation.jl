using Turing
using DifferentialEquations
using StatsPlots: plot
using LinearAlgebra
using Random
using HDF5
Random.seed!(14);

# DDE Model
function parallel_wrapper(Normal, truncated, Tsit5, solve, MethodOfSteps)
    function check_active(t::Float64, t_in_vector::Array, delay::Float64, last::Float64, is_injected::Bool)
        is_active = false
        for t_in = t_in_vector
            if (t_in + delay) < t && t < (t_in + delay + last) && is_injected
                is_active = true
            end
        end
        return is_active
    end

    function immune_response(du, u, h, p, t)
        # Free parameters 
        k6, d1, s2 = p 

        # Fixed parameters.
        t_d, t_delay, t_last, t_delay12, t_last1  = [1.556, 0.4151, 5.269, 3.98, 1.39]
        k1, k2, k3, k4, k5, _ = [0.1942, 6.04, 79.56, 1054, 5.54, 0.49]
        _, d2, d3, d4, d5, d6, d7, d8 = [11.31, 9.72, 1.25, 4.85, 0.017, 0.0381, 51.37, 0.56]
        s1, _ = [14.53, 0.3434]
        v_max=600
    
        # Current state.
        g, c, pd1, vl, vd = u
    
        # Check if treatments are active at time t
        d_cbd = check_active(t, [7], t_delay, t_last, true)
        d_12 = false
        d_cpi = false
    
        # Evaluate differential equations.
        du[1] = k1 + k2 * d_cbd - d1 * g
        du[2] = k3 + k4*h(p, t - t_d; idxs=1) - d2 * c
        du[3] = k5 - (d3+d4*g)*pd1 
        du[4] = k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
        du[5] = (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd
    
        return nothing
    end

    # Define initial-value problem
    p = (0.5, 11.0, 0.3)
    u0 = [0.0084; 9.56; 4.95; 6.65; 0]
    tspan = (0.0, 27.0)
    h(p, t; idxs::Int) = 0.0
    prob_dde = DDEProblem(immune_response, u0, h, tspan, p);

    # Solve the problem and generate a random approxiation of the solution
    selected_days = [0,7,8,9,11,14,17,20]
    num_experiments = 1
    dde_data = read_data(selected_days, num_experiments, 0.1)

    # Statistical model
    function solve_trunc(p, problem)
        predictions = solve(problem; p=p, saveat=0.1)
        pred_vol = predictions[4,:] + predictions[5,:]
        sliced_pred = pred_vol[[0, 7, 9, 11, 14, 21]*trunc(Int, 1/0.1) .+ 1]

        return sliced_pred
    end
    @model function fit_unimodal_hierarchical(data, problem, num_experiments, s, selected_days)
        if data === missing 
            data = Array{Float64}(undef, num_experiments, length(selected_days))
        end
        
        # Initialise the parameter arrays
        k6 = Vector{Float64}(undef, num_experiments)
        d1 = Vector{Float64}(undef, num_experiments)
        s2 = Vector{Float64}(undef, num_experiments)
    
        # Hyperprior distributions
        µ_k6 ~ Normal(0.5, 0.2)
        σ_k6 ~ truncated(Normal(0, 0.3); lower=0)
    
        µ_d1 ~ Normal(11, 2)
        σ_d1 ~ truncated(Normal(0, 3); lower=0)
    
        µ_s2 ~ Normal(0.3, 0.2)
        σ_s2 ~ truncated(Normal(0, 0.2); lower=0)
        
        # Regular priors
        for exp in 1:num_experiments
            k6[exp] ~ truncated(Normal(µ_k6, σ_k6); lower=0, upper=10)
            d1[exp] ~ truncated(Normal(µ_d1, σ_d1); lower=0, upper=20)
            s2[exp] ~ truncated(Normal(µ_s2, σ_s2); lower=0, upper=2)    
        end
    
        # Experimental error (σ_err)
        σ_err = 0.1 
    
        # Likelihood 
        for exp in 1:num_experiments
            p = [k6[exp], d1[exp], s2[exp]]
            sliced_pred = solve_trunc(p, problem)
    
            for i in eachindex(sliced_pred)
                data[exp, i] ~ Normal(sliced_pred[i], σ_err^2)
            end
            println("Successfully sampled from model!")
        end
    end

    # Test sampling
    chain_dde_nuts = sample(fit_unimodal_hierarchical(dde_data, prob_dde, num_experiments, 0.1, selected_days), NUTS(0.65), MCMCDistributed(),
        100, 1; progress=false)
    
    return chain_dde_nuts 
end

chain = parallel_wrapper(Normal, truncated, Tsit5, solve, MethodOfSteps)

h5open("Res/val.h5", "w") do f 
    write(f, chain)
end