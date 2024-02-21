# Functions:
# Plot more than 3 parameters

using DataFrames
using Plots
using StatsBase
using JLD2
using Statistics

plotlyjs()

# Input
abc_results = load_object("Results/abc/ABC-local-2.jld2")
true_values = [1.98814060169, 2.10643371507, 1.65262881505]
# true_values = log.([0.4037028186937715, 16.470910037502286, 0.6457205405500415])
names = ["k₆", "d₁", "s₂"]

"Plots the evoluation of the estimations for each parameter"
function plotAllEvolution()
    for i in eachindex(names)
        name = names[i]
        true_val = true_values[i]
        pops = abc_results.population[:]

        println(size.(pops))
        println(typeof(pops))
        # lower_bound = percentile.(pops,2.5)
        lower_bound = [percentile(x[:,i],2.5) for x in pops]
        upper_bound = [percentile(x[:,i],97.5) for x in pops]
        median_arr = [median(x[:,i]) for x in pops]
        mean_arr = [mean(x[:,i]) for x in pops]

        plothandle = plot()
        plot!(lower_bound; color=:red, linewidth=2)
        plot!(upper_bound; color=:red, linewidth=2)
        plot!(median_arr; color=:orange, linewidth=2)
        plot!(mean_arr; color=:purple, linewidth=2)
        xlabel!("Population Index")
        ylabel!("Paramter Value")
        title!("Approximate Posterior for $name")
        plot!(;legend=false)
        display(plothandle)
    end
end

"Plots the histogram for each parameter in the ABCResults object"
function plotAllPost()
    popn = 6
    for i in eachindex(names)
        name = names[i]
        true_val = log(true_values[i])
        final_pop = abc_results.population[8][:,i]
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

"Plots the scatter for each parameter pair"
function plotAllCorr()
    popn = 6
    combinations = [(1,2), (1,3), (2,3)]
    colours = ["#FFCCD4","#FF667D","#FF2F4E", "#D0001F", "#A20018", 
    "#ab031d","#690110"]

    for (ci, comb) in enumerate(combinations)
        i = comb[1]; j = comb[2]
        name1 = names[i]
        name2 = names[j]
        
        plothandle = plot()
        for pop_index in 1:6
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

plotAllEvolution()
plotAllPost()
# plotAllCorr()