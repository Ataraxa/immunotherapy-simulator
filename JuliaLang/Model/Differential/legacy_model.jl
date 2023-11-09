using DifferentialEquations

include("./ode_params.jl")
include("../treatments_lib.jl")

function check_active(t::Float64, t_in_vector, delay::Float64, last::Float64, is_injected::Bool)
    is_active = false
    for t_in = t_in_vector
        if (t_in + delay) < t && t < (t_in + delay + last) && is_injected
            is_active = true
        end
    end
    return is_active
end


"""
Full immune response, simulated via a delay-differential equations model
Version: 1.0 (Takuya's model)
"""
function full_immune_response(du, u, h, p, t)
    tr = p[end]
    p = p[1:21]
    
    # Model parameters.
    t_d, t_delay, t_last, t_delay12, t_last12, # 5 params
    k1, k2, k3, k4, k5, k6,
    d1, d2, d3, d4, d5, d6, d7, d8,
    s1, s2 = p
    v_max=600

    # Current state.
    g, c, pd1, vl, vd = u

    # Check if treatments are active at time t
    d_cbd = check_active(t, tr.t_in, t_delay, t_last, (tr.t_in != 0))
    d_12 = check_active(t, tr.t_in12, t_delay12, t_last12, (tr.t_in12 != 0))
    
    d_cpi = ((tr.t_inCPI) < t && (tr.t_inCPI != 0))
    # d_cpi = 0 

    # Evaluate differential equations.
    du[1] = k1 + k2 * (d_cbd + d_12) - d1 * g
    du[2] = k3 + k4*h(p, t - t_d; idxs=1) - d2 * c
    du[3] = k5 - (d3+d4*g)*pd1 
    du[4] = k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
    du[5] = (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd

    return nothing
end

struct updateParams 
    k6::Float64
    d1::Float64
    s2::Float64
end

"""
Function to perform tumour simulations in restricted parameter space.
Pass the parameters of interest in a struct, the unspecified ones will be
inferred from the average values (from GA)
"""
function restricted_simulation(
    restricted;
    tr::Treatment = CBD_IL_12_ver7,
    base::baseParams = christian,
    max_days::Float64 = 27.0
    )

    # Update the default parameters, to be passed to the DDE model
    for name in fieldnames(typeof(restricted))
        val = :($restricted.$name)
        setproperty!(base, name, eval(val))
    end
    
    # Unpack the updated structure into a vector
    vectorised_struct = Vector{Any}(undef, 22)
    for (i, name) in enumerate(fieldnames(typeof(base)))
        vectorised_struct[i] = eval(:($base.$name))
    end
    p = [vectorised_struct[1:end-1]; tr]
    u0 = [vectorised_struct[end]; 0]
    
    t_span = (0.0, max_days)
    h(p, t; idxs::Int) = 0.0

    prob_dde = DDEProblem(full_immune_response, u0, h, t_span, p)
    return prob_dde
end