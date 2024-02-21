using Distributions, DataFrames
# using Plots
using PlotlyJS: plot, scatter, attr, Layout
using LaTeXStrings

include("../../Model/treatments_lib.jl")
include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")


function do_simulations(n_iters=1_000)
    ## Variable initialisation
    all_k6  = []
    all_d1  = []
    all_s2  = []
    colours = []
    edge_len = trunc(Int64, cbrt(n_iters))

    ## Problem definition
    problem = create_problem(; 
        model="odeNfullyObs",
        max_day=27.0)
    
    params = copy(christian_true_params)
    for k6 in range(start=0.001, stop=10, length=edge_len)
        for d1 in range(start=1, stop=2, length=edge_len)
            for s2 in range(start=0.001, stop=10, length=edge_len)
                params[[11,12,21]] .= [k6, d1, s2]
                sol = solve(problem; p=params, saveat=0.1)
                # println(sol.t)
                if sol.t[end] == 27.0
                    push!(all_k6, k6)
                    push!(all_d1, d1)
                    push!(all_s2, s2)
                    push!(colours, :green)
                else
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
        font=attr( 
            size =15
        ),

        scene=attr(
            # xaxis_type="log",
            # yaxis_type="log",
            # zaxis_type="log",

            xaxis_title="k6",
            yaxis_title="d1",
            zaxis_title="s2",
            
            # xaxis_tickvals=[0.1, 1,10,100],
            # yaxis_tickvals=[0.1, 1,10,100],
            # zaxis_tickvals=[0.01,0.1, 1,10,100],
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

do_simulations(1000)
