using GLMakie
using Match
import Makie.plot

# Fetch :Chains object from file
chain = h5open("Results/hpc-individual-1-1.h5", "r") do f
    read(f, Chains)
end

# Declare true values for each parameter, to draw vlines
true_p = Dict(
    :ln_k6 => -0.6171,
    :ln_d1 => 2.3201,
    :ln_s2 => -0.8896
)

function Makie.plot(ch::Chains)
    fig = Figure()
    for (ind, param) in enumerate(ch.name_map.parameters)
        global ax = Axis(fig[ind, 1], title=string(param))
        for (ind2, datavec) in enumerate(eachcol(getindex(ch, param).data))
            if ind2 in [4] # Ignore badly mixed chains 
                continue
            end
            # Get current default colorpalette
            colors = Makie.current_default_theme().attributes[:palette][][:color][]
            density!(ax, datavec; color=(:black, 0.0),
                strokearound=true,
                strokewidth=2,
                strokecolor=colors[ind2%length(colors)],
            ) 
        end
        # Add vertical line based on true value 
        vlines!(true_p[param]; color=:red)

        # Change xlim bounds if necessary (manually)
        @match param begin
            :ln_s2 => xlims!.((ax), -5, 0)
            _ => ""
        end
        # xlims!.((ax), max(-2,xmin), min(2, xmax))
    end
    display(fig)
    return fig
end

Makie.plot(chain)
