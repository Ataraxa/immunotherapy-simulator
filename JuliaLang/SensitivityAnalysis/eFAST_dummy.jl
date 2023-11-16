using GlobalSensitivity
using DifferentialEquations
using Distributions
using Plots
using StatsPlots: groupedbar

# Settings
np = 3 # number of parameters to be analysed

# Wrap the problem inside a function that solves it with a given set of params 
wrapped_response = function (p)
    return p[1]*10+p[2]
end 

# Lower and upper bounds of the solution 
vector = [1, 1, 1]
width_factor = 0.5
lb = (1-width_factor)*vector
ub = (1+width_factor)*vector

# Labels must be in same order as parameter vector
labels = [
    "a", "b", "c"
    ]
res_sens = gsa(wrapped_response, eFAST(),[[lb[i], ub[i]] for i=1:np], samples=500)
# print(size(res_sens.S1))
# scatter(1:np, (res_sens.S1+res_sens.ST)'; ticks=(1:np,labels))
println("n₁ : S1=$(res_sens.S1[1, 3]) | ST=$(res_sens.ST[1, 3])")
groupedbar([res_sens.S1; res_sens.ST]'; 
    bar_position=:stack,
    xticks=(1:np,labels),
    group=repeat(["First-Order", "Total-Order"], inner=np))

