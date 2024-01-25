using DifferentialEquations
using Match

include("./treatments_lib.jl")
include("./Bayesian/priors.jl")

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
Function that returns the default, fully functional DDE Problem
Inputs:
    - model: the tumour DDE model used to perform simulations
    - treatment: the treatment to use for the simulation

Outputs:
    - ::DDEProblem, to be solved later
"""
function model_factory(; model="takuya", treatment::Treatment=CBD_IL_12_ver7)
    immune_resp::Function = (x->print(x))
    
    @match model begin
       "takuya" => begin 
            immune_resp = function(du, u, h, p, t)
                tr = treatment
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
            
                return du
            end
        end

        "w/feedback" => begin 
        immune_resp = function(du, u, h, p, t)
            tr = treatment
            p = p[1:21]
            
            # Model parameters.
            t_d, t_delay, t_last, t_delay12, t_last12, # 5 params
            k1, k2, k3, k4, k5, k6,
            d1, d2, d3, d4, d5, d6, d7, d8,
            s1, s2 = p
            v_max=600

            new = 8
        
            # Current state.
            g, c, pd1, vl, vd = u
        
            # Check if treatments are active at time t
            d_cbd = check_active(t, tr.t_in, t_delay, t_last, (tr.t_in != 0))
            d_12 = check_active(t, tr.t_in12, t_delay12, t_last12, (tr.t_in12 != 0))
            
            d_cpi = ((tr.t_inCPI) < t && (tr.t_inCPI != 0))
            # d_cpi = 0 
        
            # Evaluate differential equations.
            du[1] = k1 + k2 * (d_cbd + d_12) - d1 * g + new * g
            du[2] = k3 + k4*h(p, t - t_d; idxs=1) - d2 * c
            du[3] = k5 - (d3+d4*g)*pd1 
            du[4] = k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
            du[5] = (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd
        
            return du
        end
    end
    end

    return immune_resp
end

function create_problem(; 
        model="takuya", 
        treatment::Treatment=CBD_IL_12_ver7, 
        params=christian_true_params,
        max_day::Float64 = 27.0
        )
    
    model=model_factory(model=model, treatment=treatment)
    p = params[1:21]
    u0 = [params[22:end]; 0]
    h(p, t; idxs::Int) = 0.0
    t_span = (0.0, max_day)

    return DDEProblem(model, u0, h, t_span, p)
end