using Distributions, DataFrames
# using Plots
using PlotlyJS: plot, scatter, attr, Layout
using LaTeXStrings

include("../../Model/treatments_lib.jl")
include("../../Model/Differential/ode_core.jl")
include("../../Model/Differential/ode_params.jl")
include("../../Model/Differential/ode_restricted.jl")


function do_simulations(n_iters=1_000, θ=25)
    ## Variable initialisation
    all_k6  = []
    all_d1  = []
    all_s2  = []
    colours = []
    edge_len = trunc(Int64, cbrt(n_iters))

    ## Problem definition
    problem = create_problem(; 
        model="takuya",
        param_struct=christian,
        max_day=100.0)
    
    for k6 in exp10.(range(-1, stop=1, length=edge_len))
        for d1 in exp10.(range(0, stop=2, length=edge_len))
            for s2 in exp10.(range(-2, stop=0, length=edge_len))
                p, _ = repack_params(updateParams3(k6, d1, s2), christian)
                sol = solve(problem; p=p, saveat=0.1)

                total_tumour = sol[4,:] + sol[5,:]

                if mean(total_tumour[end-20:end]) < θ
                    push!(all_k6, k6)
                    push!(all_d1, d1)
                    push!(all_s2, s2)
                    push!(colours, :red)
                else
                    push!(all_k6, k6)
                    push!(all_d1, d1)
                    push!(all_s2, s2)
                    push!(colours, :green)
                end
            end
        end
    end

    println(typeof(all_k6))

    layout = Layout(
        font=attr( 
            size =15
        ),

        scene=attr(
            xaxis_type="log",
            yaxis_type="log",
            zaxis_type="log",

            xaxis_title="k6",
            yaxis_title="d1",
            zaxis_title="s2",
            
            xaxis_tickvals=[0.1, 1,10,100],
            yaxis_tickvals=[0.1, 1,10,100],
            zaxis_tickvals=[0.01,0.1, 1,10,100],
            ),
        
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