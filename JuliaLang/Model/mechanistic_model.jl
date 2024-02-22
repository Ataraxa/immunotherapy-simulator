using DifferentialEquations
using Match
using ForwardDiff

include("./treatments_lib.jl")
include("./Bayesian/priors.jl")

function check_active(t::Union{Float64, ForwardDiff.Dual}, t_in_vector, delay::Float64, last::Float64, is_injected::Bool)
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
            end # closes model
        end # closes begin loop

        "w/feedback" => begin 
            immune_resp = function(du, u, h, p, t)
                tr = treatment
                
                # Model parameters.
                t_d, t_delay, t_last, t_delay12, t_last12, # 5 params
                k1, k2, k3, k4, k5, k6,
                d1, d2, d3, d4, d5, d6, d7, d8,
                s1, s2,
                n1 = p
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
                du[1] = k1 + k2 * (d_cbd + d_12) - d1 * g + n1 * g
                du[2] = k3 + k4*h(p, t - t_d; idxs=1) - d2 * c
                du[3] = k5 - (d3+d4*g)*pd1 
                du[4] = k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
                du[5] = (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd
            
                return du
            end # closes function
        end #closes begin loop

        "odeNnon" => begin
            immune_resp = function(du, u, p, t)
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
            
                # Evaluate differential equations.
                du[1] = k1 + k2 * (d_cbd + d_12) - d1 * g
                du[2] = k3 + k4*g - d2 * c
                du[3] = k5 - (d3+d4*g)*pd1 
                du[4] = k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
                du[5] = (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd
            
                return du
            end # closes function block
        end # closes specific case statement

        "odeNfullyObs" => begin
            immune_resp = function(du, u, p, t)
                tr = treatment
                p = p[1:21]
                
                # Model parameters.
                t_d, t_delay, t_last, t_delay12, t_last12, # 5 params
                k1, k2, k3, k4, k5, k6,
                d1, d2, d3, d4, d5, d6, d7, d8,
                s1, s2 = p
                v_max=600
            
                # Current state.
                g, c, pd1, vl = u
            
                # Check if treatments are active at time t
                d_cbd = check_active(t, tr.t_in, t_delay, t_last, (tr.t_in != 0))
                d_12 = check_active(t, tr.t_in12, t_delay12, t_last12, (tr.t_in12 != 0))
                
                d_cpi = ((tr.t_inCPI) < t && (tr.t_inCPI != 0))
            
                # Evaluate differential equations.
                du[1] = k1 + k2 * (d_cbd + d_12) - d1 * g
                du[2] = k3 + k4*pd1 - d2 * c
                du[3] = k5 - (d3+d4*g)*pd1 
                du[4] = k6*(1-(vl)/v_max)*vl - (d5 + d7*g + d6*c - s1*pd1)
                return du
            end # closes function block
        end # closes specific case statement

        "ddeNfullyObs" => begin
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
                g, c, pd1, vl = u
            
                # Check if treatments are active at time t
                d_cbd = check_active(t, tr.t_in, t_delay, t_last, (tr.t_in != 0))
                d_12 = check_active(t, tr.t_in12, t_delay12, t_last12, (tr.t_in12 != 0))
                
                d_cpi = ((tr.t_inCPI) < t && (tr.t_inCPI != 0))
            
                # Evaluate differential equations.
                du[1] = k1 + k2 * (d_cbd + d_12) - d1 * g
                du[2] = k3 + k4*h(p, t - t_d; idxs=1) - d2 * c
                du[3] = k5 - (d3+d4*g)*pd1 
                du[4] = k6*(1-(vl)/v_max)*vl - (d5 + d7*g + d6*c - s1*pd1)
                return du
            end # closes function block
        end # closes specific case statement

        "predatorPrey" => begin
            immune_resp = function(du, u, p, t)
                # Model parameters.
                α, β, γ, δ = p
                # Current state.
                x, y = u
            
                # Evaluate differential equations.
                du[1] = (α - β * y) * x # prey
                du[2] = (δ * x - γ) * y # predator
                return du
            end
        end # closes specific case statement
    end #closes match statement

    return immune_resp
end # closes factory function

function create_problem(; 
        model="takuya", 
        treatment::Treatment=CBD_IL_12_ver7, 
        params=christian_true_params,
        max_day::Float64 = 27.0
        )
    
    model_obj=model_factory(model=model, treatment=treatment)
    p = params[1:21]
    u0 = [params[22:end]; 0]
    h(p, t; idxs::Int) = 0.0
    t_span = (0.0, max_day)

    if model=="odeNnon"
        return ODEProblem(model_obj, u0, t_span, p)
    elseif model=="odeNfullyObs"
        return ODEProblem(model_obj, u0[1:4], t_span, p) # no vd
    elseif model=="ddeNfullyObs"
        return DDEProblem(model_obj, u0[1:4], h, t_span, p)
    elseif model=="predatorPrey"
        u0 = [1.0, 1.0]
        p = [1.5, 1.0, 3.0, 1.0]
        tspan = (0.0, 10.0)
        return ODEProblem(model_obj, u0, tspan, p)
    else
        return DDEProblem(model_obj, u0, h, t_span, p)
    end
end
