# Functions:
# Plot more than 3 parameters

using DataFrames
using Plots
using StatsBase
using JLD2
# Plots.scalefontsizes(1/1.5)

# Input
abc_results = load_object("Results/abc/ABC-local-0.jld2")
true_values = log.([1.98814060169, 2.10643371507, 1.65262881505])
names = ["k₆", "d₁", "s₂"]


function plotAllEvolution()
    for i in eachindex(names)
        name = names[i]
        true_val = true_values[i]
        final_pop = abc_results.population[6][:,i]
        distance = abc_results.distances[6]
        threshold = abc_results.threshold_schedule[6]

        lower_bound = percentile(final_pop,2.5)
        upper_bound = percentile(final_pop,97.5)

        global df = DataFrame([final_pop, distance], ["p", "dist"])
        plothandle = plot()
        histogram!(df.p, bins=range(minimum(df.p),maximum(df.p), length=50))
        vline!([true_val], linewidth=5)
        vline!([lower_bound upper_bound]; fill=true,color=:red, linewidth=5)
        xlabel!("Parameter Value")
        ylabel!("Number of Accepted Particles")
        title!("Approximate Posterior for $name")
        plot!(;legend=false)
        display(plothandle)
    end
end

function plotAllPost()
    popn = 7
    for i in eachindex(names)
        name = names[i]
        true_val = true_values[i]
        final_pop = abc_results.population[popn][:,i]
        distance = abc_results.distances[popn]
        threshold = abc_results.threshold_schedule[popn]

        lower_bound = percentile(final_pop,2.5)
        upper_bound = percentile(final_pop,97.5)

        global df = DataFrame([final_pop, distance], ["p", "dist"])
        plothandle = plot()
        histogram!(df.p, bins=range(minimum(df.p),maximum(df.p), length=50))
        vline!([true_val], linewidth=5)
        vline!([lower_bound upper_bound]; fill=true,color=:red, linewidth=5)
        xlabel!("Parameter Value")
        ylabel!("Number of Accepted Particles")
        title!("Approximate Posterior for ln($name)")
        xlims!((-1., 2.5))
        plot!(;legend=false)
        display(plothandle)
    end
end

function plotAllCorr()
    popn = 7
    combinations = [(1,2), (1,3), (2,3)]
    colours = ["#FFCCD4","#FF667D","#FF2F4E", "#D0001F", "#A20018", 
    "#ab031d","#690110"]

    for (ci, comb) in enumerate(combinations)
        i = comb[1]; j = comb[2]
        name1 = names[i]
        name2 = names[j]
        
        plothandle = plot()
        for pop_index in 1:7
            pop1 = abc_results.population[pop_index][:,i]
            pop2 = abc_results.population[pop_index][:,j]
        scatter!(pop1, pop2; color=colours[pop_index], markersize=5)
        end
        xlabel!("Parameter Value of ln($name1)")
        ylabel!("Parameter Value of ln($name2)")
        # title!("Correlation")
        plot!(;legend=false)
        display(plothandle)
    end
end
plotAllPost()