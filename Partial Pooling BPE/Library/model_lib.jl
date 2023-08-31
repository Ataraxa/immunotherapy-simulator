# File containing the necessary objects to simulate the immune response.

"""
Short function to check is a treament (characterised by the params)is active at
time t.
"""
function check_active(t::Float64, t_in_vector::Array, delay::Float64, last::Float64, is_injected::Bool)
    is_active = false
    for t_in = t_in_vector
        if (t_in + delay) < t && t < (t_in + delay + last) && is_injected
            is_active = true
        end
    end
    return is_active
end


"""
DDE model as a ODEProblem object.
"""
function immune_response(du, u, h, p, t)
    # Model parameters.
    t_d, t_delay, t_last, t_delay12, t_last1  = p[1]
    k1, k2, k3, k4, k5, k6 = p[2]
    d1, d2, d3, d4, d5, d6, d7, d8 = p[3]
    s1, s2 = p[4]
    v_max=600

    # Current state.
    # print(u)
    g, c, p, vl, vd = u

    # Check if treatments are active at time t
    d_cbd = check_active(t, [7], t_delay, t_last, true)
    d_12 = false
    d_cpi = false

    # Evaluate differential equations.
    du[1] = k1 + k2 * d_cbd - d1 * g
    du[2] = k3 + k4*h(p, t - t_d; idxs=1) - d2 * c
    du[3] = k5 - (d3+d4*g)*p 
    du[4] = k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*p*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
    du[5] = (d5 + (d6*c/(1+s1*p*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd

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
    p=[pt, pk, pd, ps]

    return u0, p
end