using Plots

include("./continuation_lib.jl")

function simplified_tumour(u, p)
    vₘₐₓ = 600
    v = u[1]
    dv = (1-v/vₘₐₓ)*v - p*v
    return [dv]
end

function simplified_jacob(u, p)
    vₘₐₓ = 600
    v = u[1]
    dv = (1-p) - 2*v / vₘₐₓ
    return [dv]
end

pmin = 0.00  ; pmax = 100; δ = 0.9
x0 = [10.0]  ; p0 = 0
dx0 = [0.1]; dp0 = 0.01

xs, ps, stability = continuation(simplified_tumour, simplified_jacob, x0, p0;
    pmin, pmax, dp0, dx0
)

colors = [s ? :blue : :red for s in stability]
p = scatter(ps, [x[1] for x in xs]; color = colors, markerstrokecolor = colors, 
    xlabel = "p", ylabel = "x", label = "")