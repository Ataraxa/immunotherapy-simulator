using Distributions, DataFrames
using PlotlyJS: plot, scatter, attr, Layout

include("../../Model/Differential/ode_restricted.jl")
include("../../Model/treatments_lib.jl")

function do_simulations(n_iters=1_000, θ=25)
    all_k6 = Vector{Float64}(undef, n_iters)
    all_d1 = Vector{Float64}(undef, n_iters)
    all_s2 = Vector{Float64}(undef, n_iters)
    colours = Vector{Symbol}(undef, n_iters)

    for i in 1:n_iters
        k6 = rand(Uniform(0.0, 1.0 ))
        d1 = rand(Uniform(7.0, 15.0))
        s2 = rand(Uniform(0.0, 1.0 ))
        
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
    matrix = [reshape(all_k6, 1, :); reshape(all_d1, 1, :); reshape(all_s2, 1, :)]
    return plot(scatter(
        x=all_k6,
        y=all_d1,
        z=all_s2,
        mode="markers",
        marker=attr(
            size=12,
            color=colours,                # set color to an array/list of desired values
            colorscale="Viridis",   # choose a colorscale
            opacity=0.8
        ),
        type="scatter3d"
    ), Layout(margin=attr(l=0, r=0, b=0, t=0)))
end

do_simulations(1000, 50)