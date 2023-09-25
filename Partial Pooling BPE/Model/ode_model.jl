# File containing the necessary objects to simulate the immune response.

# """
# Short function to check is a treament (characterised by the params)is active at
# time t.
# """

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
    t_d, t_delay, t_last, t_delay12, t_last1, # 5 params
    k1, k2, k3, k4, k5, k6,
    d1, d2, d3, d4, d5, d6, d7, d8,
    s1, s2 = p
    v_max=600

    # Current state.
    g, c, pd1, vl, vd = u

    # Check if treatments are active at time t
    d_cbd = check_active(t, tr["t_in"], t_delay, t_last, tr["active_cbd"])
    d_12 = check_active(t, tr["t_in12"], t_delay, t_last, tr["active_il12"])
    
    # d_cpi = ((tr["t_inCPI"]) < t && tr["active_cpi"])
    d_cpi = 0 

    # Evaluate differential equations.
    du[1] = k1 + k2 * (d_cbd + d_12) - d1 * g
    du[2] = k3 + k4*h(p, t - t_d; idxs=1) - d2 * c
    du[3] = k5 - (d3+d4*g)*pd1 
    du[4] = k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
    du[5] = (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd

    return nothing
end

"""
Full immune response, simulated via a deay-differential equations model
Version: 1.1 (Adding the positive feedback loop)
"""
function full_with_feedback(du, u, h, p, t)
    # Model parameters.
    t_d, t_delay, t_last, t_delay12, t_last1,
    k1, k2, k3, k4, k5, k6,
    d1, d2, d3, d4, d5, d6, d7, d8,
    s1, s2,
    b1 = p
    v_max=600

    # Current state.
    g, c, pd1, vl, vd = u

    # Check if treatments are active at time t
    d_cbd = check_active(t, [7], t_delay, t_last, true)
    d_12 = false
    d_cpi = false

    # Evaluate differential equations.
    du[1] = k1 + k2 * d_cbd + (b1 - d1) * g
    du[2] = k3 + k4*h(p, t - t_d; idxs=1) - d2 * c
    du[3] = k5 - (d3+d4*g)*pd1 
    du[4] = k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
    du[5] = (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd

    return nothing
end

"""Default values getter"""
function get_default_values()
    # Define initial-value problem.
    u0 = [0.0084; 9.56; 4.95; 6.65; 0]

    pt = [1.556, 0.4151, 5.269, 3.98, 1.39]
    pk = [0.1942, 6.04, 79.56, 1054, 5.54, 0.49]
    pd = [11.31, 9.72, 1.25, 4.85, 0.017, 0.0381, 51.37, 0.56]
    ps = [14.53, 0.3434]
    p=[pt; pk; pd; ps] # Combine all in a single vector

    return u0, p
end