using BifurcationKit

function simplified_tumour(dv, v, p, t=0)
    vₘₐₓ = 600
    dv = v*(1 - p - v/vₘₐₓ)
    dv
end

parm = 0.5
v0 = 6.0

prob = BifurcationProblem(simplified_tumour, v0, parm)

opts_br = ContinuationPar(p_min = 0.0, p_max = 10.0, ds = 0.04, dsmax = 0.05,)

br = continuation(prob, PALC(tangent=Bordered()), opts_br; normC = norminf)

scene = plot(br, plotfold=false, markersize=3, legend=:topleft)