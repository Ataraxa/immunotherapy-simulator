#=
Julai implementation the Approximate Bayesian Compuation (ABC) tutorial from the
paper https://doi.org/10.1016/j.jmp.2012.02.005
=#
using Distributions
using ApproxBayes
using Plots: plot, xlims!, Axis
using StatsPlots: plot

n = 100
function dist(params, constants,  targetdata)
    p = params[1]
    simulation = [(rand(Uniform(0,1)) >= p) for _ in 1:n]
    Y = sum(simulation)
    X = sum(targetdata)
    distance = 1/n * abs(Y - X), 1
    return distance
end

targetdata = [(rand() >= 0.7) for _ in 1:n]
Y=sum(targetdata)

setup = ABCRejection(
    dist,
    1,
    1e-10,
    Prior([Beta(1,1)])
)

smc = runabc(setup, targetdata)
a = plot(smc; title = "géfékk")
a = plot!(Beta(n-Y + 1, Y + 1))


