using DifferentialEquations

"""
Adapts the base DDE model for Bayesian fitting by reducing the parameter space.
Parameters of interest (given by sensitivity analysis) are free parameters,
the remaining ones are set to the average value determined by GA optimisation
"""

# Import the full model
include("ode_model.jl")

# Adapts the model and packs it into a DDE problem
function adapt_dde_space()
    u0, _ = get_default_values()
    p = [0.6, 11, 0.4]
    t_span = (0.0, 27.0)
    h(p, t; idxs::Int) = 0.0
    
    prob_dde = DDEProblem(bayesian_immune_response, u0, h, t_span, p)

    return prob_dde
end

function bayesian_immune_response(du, u, h, p, t)
    # Free parameters 
    k6, d1, s2 = p 

    # Fixed parameters
    _, def_params = get_default_values()

    t_d, t_delay, t_last, t_delay12, t_last1  = def_params[1]
    k1, k2, k3, k4, k5, _ = def_params[2]
    _, d2, d3, d4, d5, d6, d7, d8 = def_params[3]
    s1, _ = def_params[4]
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
end