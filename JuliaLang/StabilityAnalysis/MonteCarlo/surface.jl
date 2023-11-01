using Distributions, DataFrames
using PlotlyJS: plot, scatter, attr, Layout
using LaTeXStrings

include("../../Model/Differential/ode_restricted.jl")
include("../../Model/treatments_lib.jl")

function do_simulations(n_iters=1_000, θ=25)
    all_k6  = Vector{Float64}(undef, n_iters)
    all_d1  = Vector{Float64}(undef, n_iters)
    all_s2  = Vector{Float64}(undef, n_iters)
    colours = Vector{Symbol}(undef, n_iters)
    
    i = 0
    for k6 in exp10.(range(-1, stop=1, length=10))
        for d1 in exp10.(range(0, stop=2, length=10))
            for s2 in exp10.(range(-2, stop=0, length=10))
                i += 1
                all_k6[i] = k6
                all_d1[i] = d1
                all_s2[i] = s2

                update = updateParams(k6, d1, s2)
                prob = restricted_simulation(update)
                sol = solve(prob; saveat=0.1)

                total_tumour = sol[4,:] + sol[5,:]
                is_responder = (total_tumour[end] < θ) ? :green : :red
                colours[i] = is_responder
            end
        end
    end

    layout = Layout(
        scene=attr(
            xaxis_type="log",
            yaxis_type="log",
            zaxis_type="log",

            xaxis_title="k6",
            yaxis_title="d1",
            zaxis_title="s2",
            )
        )

    return plot(scatter(
        x=all_k6,
        y=all_d1,
        z=all_s2,
        mode="markers",
        marker=attr(
            size=20,
            color=colours,                # set color to an array/list of desired values
            colorscale="Viridis",   # choose a colorscale
            opacity= (colours == :red) ? 0 : 1
        ),
        type="scatter3d"
    ), layout)
end

do_simulations(1000, 50)