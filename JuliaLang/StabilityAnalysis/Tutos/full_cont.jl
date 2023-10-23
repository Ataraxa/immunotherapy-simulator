using Revise, Parameters, Plots
using BifurcationKit
const BK = BifurcationKit
using Distributions: Uniform

function simplified_tumour(du, u, p, t = 0)
    # vₘₐₓ = 600
    v = u[1]
    # du[1] = v*(1 - p - v/vₘₐₓ)
    # du[1] = p + v^2 # Saddle node 
    du[1] = v*(p-v) # Transcritical
    return du 
end

parm = -4.0
v0 = [rand(Uniform(-4.0, -2.0))]
# v0 = -4.0

prob = BifurcationProblem(simplified_tumour, v0, parm,
    record_from_solution = (x, p) -> (v = x[1]))

opts_br = ContinuationPar(p_min = -5.0, p_max = 5.0, 
    # To ensure smooth continuation
    ds = 0.04, dsmax = 0.1,)

br = continuation(prob, PALC(tangent=Secant()), opts_br; normC=norminf)

scene = plot(br, plotfold=true, markersize=3, legend=:topleft)