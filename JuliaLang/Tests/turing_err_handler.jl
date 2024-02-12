@model function test_model(
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

    if pred.t[end] != 27.0
        println.(float_p)
        println("___________________________")
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