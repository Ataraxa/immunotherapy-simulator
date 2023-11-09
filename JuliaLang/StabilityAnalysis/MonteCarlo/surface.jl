using Distributions, DataFrames
# using Plots
using PlotlyJS: plot, scatter, attr, Layout
using LaTeXStrings

include("../../Model/Differential/legacy_model.jl")
include("../../Model/treatments_lib.jl")

function do_simulations(n_iters=1_000, θ=25)
    all_k6  = []
    all_d1  = []
    all_s2  = []
    colours = []
    
    for k6 in exp10.(range(-1, stop=1, length=10))
        for d1 in exp10.(range(0, stop=2, length=10))
            for s2 in exp10.(range(-2, stop=0, length=10))
                update = updateParams(k6, d1, s2)
                prob = restricted_simulation(update; max_days=100.0)
                sol = solve(prob; saveat=0.1)

                total_tumour = sol[4,:] + sol[5,:]

                if mean(total_tumour[end-20:end]) < θ
                    push!(all_k6, k6)
                    push!(all_d1, d1)
                    push!(all_s2, s2)
                    push!(colours, :red)
                end
            end
        end
    end

    println(typeof(all_k6))

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
            size=10,
            color=colours,        # set color to an array/list of desired values
            colorscale="Viridis",
            opacity=0.8 # choose a colorscale
        ),
        type="scatter3d"
    ), layout)
end

do_simulations(1000, 10)