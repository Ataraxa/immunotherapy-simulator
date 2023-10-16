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
    t_d, t_delay, t_last, t_delay12, t_last12, # 5 params
    k1, k2, k3, k4, k5, k6,
    d1, d2, d3, d4, d5, d6, d7, d8,
    s1, s2 = p
    v_max=600

    # Current state.
    g, c, pd1, vl, vd = u

    # Check if treatments are active at time t
    d_cbd = check_active(t, tr["t_in"], t_delay, t_last, tr["active_cbd"])
    d_12 = check_active(t, tr["t_in12"], t_delay12, t_last12, tr["active_il12"])
    
    d_cpi = ((tr["t_inCPI"]) < t && tr["active_cpi"])
    # d_cpi = 0 

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

function get_christian()
    p=[1.4700, 0.3800, 5.4100, 3.7400, 1.3600, 0.1700, 5.3000, 83.9000,
     1292.2, 4.8500, 0.4900, 9.9700, 11.0500, 1.0900,
     5.8000, 0.0200, 0.0400, 51.9400, 0.6100, 16, 0.3500]

    u0=[0.0090, 10, 4.4000, 7.0600, 0]

    return u0, p
end

function get_temp_test()
    p = [1.5560, 0.4151, 5.2695, 3.9868, 1.3903, 0.1942, 6.0401, 79.5663, 
    1.0546e+03, 5.4378, 0.4901, 11.3173, 9.7233, 1.2598, 4.8524, 0.0171, 
    0.0381, 51.3708, 0.5607, 14.5359, 0.3434] 
    
    u0 = [0.0084, 9.5608, 4.9575, 6.6580, 0]
    
    return u0, p
end

function get_temp_test2()

    p = [4.9487, 4.2399, 9.9006, 0.4167, 8.3883, 3.0158, 0.9189, 154.1321,
    874.7748, 0.4008, 3.4200, 6.4546, 13.7066, 9.5774, 5.6896, 2.1396, 7.7700,
    84.1471, 0.4819, 48.1169, 8.5841]

    u0 = [0.7994, 9.7683, 1.6256, 5.6674, 0]
    
    return u0, p
end