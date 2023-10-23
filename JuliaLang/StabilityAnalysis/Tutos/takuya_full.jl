using Revise, Parameters, Plots
using BifurcationKit
const BK = BifurcationKit

include("../Model/Differential/ode_params.jl")

# vector field
function TMvf!(dz, z, p, t = 0)
	@unpack t_d, t_delay, t_last, t_delay12, t_last12, 
    k1, k2, k3, k4, k5, k6, 
    d1, d2, d3, d4, d5, d6, d7, d8, 
    s1, s2 = p
	E, x, u = z
	SS0 = J * u * x * E + E0
	SS1 = α * log(1 + exp(SS0 / α))
	dz[1] = (-E + SS1) / τ
	dz[2] =	(1.0 - x) / τD - u * x * E
	dz[3] = (U0 - u) / τF +  U0 * (1.0 - u) * E
	dz
end

# parameter values
par_tm = christian

# initial condition
z0 = [0.238616, 0.982747, 0.367876]

# Bifurcation Problem
prob = BifurcationProblem(TMvf!, z0, par_tm, (@lens _.k6);
	record_from_solution = (x, p) -> (E = x[1], x = x[2], u = x[3]),)

opts_br = ContinuationPar(p_min = -10.0, p_max = -0.9,
    # parameters to have a smooth continuation curve
    ds = 0.04, dsmax = 0.05,)

# continuation of equilibria
br = continuation(prob, PALC(tangent=Bordered()), opts_br; normC = norminf)
    
scene = plot(br, plotfold=true, markersize=3, legend=:topleft)