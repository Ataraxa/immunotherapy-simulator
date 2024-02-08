using Plots
using StatsBase

abc_results = load_object("Results/abc/ABC-local-0.jld2")

for i in 1:3
    a = abc_results.population[1:end-1]
    x1=percentile.([f[:,i] for f in a], 2.5)
    x2=percentile.([f[:,i] for f in a], 97.5)
    x3=median.([f[:,i] for f in a])

    p1 = plot()
    plot!(1:7,x1, color=:red)
    plot!(1:7,x2, color=:red)
    plot!(1:7,x3, color=:blue)
    display(p1)
end
