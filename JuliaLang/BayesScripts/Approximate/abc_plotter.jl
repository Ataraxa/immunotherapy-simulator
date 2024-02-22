# Functions:
# Plot more than 3 parameters

using DataFrames
using Plots
using StatsBase
using JLD2
using Statistics

plotlyjs()

# Input
abc_results = load_object("Results/abc/ABC-hpc-1.jld2")
# true_values = [1.98814060169, 2.10643371507, 1.65262881505]
# names = ["k₆", "d₁", "s₂"]

# For large inference
true_values = [1.872824160505723, # t_d 
0.487923336960266, # t_delay 
4.893475496158970, # t_last 
3.697285769720824, # t_delay1 
1.075411285414352, # t_last12
# ⬇ index 6
0.221936172352179,  #k1
6.081963924603686,  #k2
74.60149998881774,  #k3
929.5659753695654,  #k4
5.815052050675676,  #k5
0.539494697127057,  #k6
# ⬇ index 12
10.177238034144443, # d1 
10.610875554141064, # d2 
1.3262458309759877, # d3
4.4860998573217170, # d4
0.0160245808033653, # d5
0.0341621446708329, # d6
59.611074450040185, # d7
0.5697851077620799, # d8
# ⬇ index 20
14.068450976326906, # s1
0.4108126097661595, # s2
# ⬇ index 22
0.007911968983770, # g (IFNγ)
8.851474387488993, # c (CD8+)
5.999229998549251, # p (PD-1)
5.573047884761726  # v (Living tumour)
]
names = ["t_d", "t_delay", "t_last", "t_delay12", "t_last12", # 5 params
"k1", "k2", "k3", "k4", "k5", "k6",
"d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8",
"s1", "s2" ]

"Plots the evoluation of the estimations for each parameter"
function plotAllEvolution()
    for i in eachindex(names)
        name = names[i]
        # true_val = true_values[i]
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