using Revise, DDEBifurcationKit, Parameters, LinearAlgebra, Plots
using BifurcationKit
const BK = BifurcationKit

include("../Model/treatments_lib.jl")
include("../Model/Differential/ode_params.jl")

function immune_resp(u, ud, p::baseParams)
   tr::Treatment = CBD_IL_12_ver7
   
   # Model parameters.
   @unpack t_d, t_delay, t_last, t_delay12, t_last12, # 5 params
   k1, k2, k3, k4, k5, k6,
   d1, d2, d3, d4, d5, d6, d7, d8,
   s1, s2 = p

   v_max=600

   # Current state.
   g, c, pd1, vl, vd = u
   gd, cd, pd1d, vld, vdd = ud[1] # most of the delayed values are useless

   # Check if treatments are active at time t
   d_cbd = check_active(t, tr.t_in, t_delay, t_last, (tr.t_in != 0))
   d_12 = check_active(t, tr.t_in12, t_delay12, t_last12, (tr.t_in12 != 0))
   d_cpi = ((tr.t_inCPI) < t && (tr.t_inCPI != 0))

   # Evaluate differential equations.
   [
      k1 + k2 * (d_cbd + d_12) - d1 * g
      k3 + k4*gd - d2 * c
      k5 - (d3+d4*g)*pd1 
      k6*(1-(vl+vd)/v_max)*vl - (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl
      (d5 + (d6*c/(1+s1*pd1*(1-d_cpi)) + d7*g)/(1+s2*(vl+vd)))*vl - d8*vd 
   ]
end

delaysF(par) = [par.t_d]

pars = christian
x0 = [pars.u0; 0]

prob = ConstantDDEBifProblem(immune_resp, delaysF, x0, pars, (@lens _.s2))

optn = NewtonPar(verbose = true, eigsolver = DDE_DefaultEig())
opts = ContinuationPar(p_max = 13., p_min = 0., newton_options = optn, ds = -0.01, detect_bifurcation = 3, nev = 5, dsmax = 0.2, n_inversion = 4)
br = continuation(prob, PALC(), opts; verbosity = 1, plot = true, bothside = true, normC = norminf)

scene = plot(br)