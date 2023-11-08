using GlobalSensitivity
using DifferentialEquations
using Distributions
using Plots

include("../Model/Differential/ode_core.jl")
include("../Model/Differential/ode_params.jl")
include("../CommonLibrary/struct_manipulation.jl")

# Settings
np = 23 # number of parameters to be analysed

# Define the problem object
t_span = (0.0, 25.0)
h(p, t; idxs::Int) = 0.0
p, u0 = struct_split(christian)
problem = DDEProblem(full_with_feedback, u0, h, t_span, p)

# Wrap the problem inside a function that solves it with a given set of params 
wrapped_response = function (params_vector)
    # Assume vector = [p; vl0]
    p = params_vector[1:end-1] 
    u0 = [0.0084; 9.56; 4.95; params_vector[end]; 0]
    prob1 = remake(problem;p=p, u0=u0)
    sol = solve(prob1;saveat=0.1)
    return(sol[4,end] + sol[5,end])
end 

# Lower and upper bounds of the solution 
vector = [p; 2.0; u0[4];]
width_factor = 10 
lb = (1/width_factor)*vector
ub = width_factor*vector

# Labels must be in same order as parameter vector
labels = [
    "tem", "td", "tl", "NA", "NA",
    "k1", "k2", "k3", "k4", "k5", "k6",
    "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8",
    "s1", "s2", "v0"
    ]
res_sens = gsa(wrapped_response, eFAST(),[[lb[i], ub[i]] for i=1:np], samples=500)
scatter(1:np, (res_sens.S1+res_sens.ST)'; ticks=(1:np,labels))
