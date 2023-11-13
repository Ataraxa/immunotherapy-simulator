using GlobalSensitivity
using DifferentialEquations
using Distributions
using Plots
using StatsPlots: groupedbar

include("../Model/treatments_lib.jl")
include("../Model/Differential/ode_core.jl")
include("../Model/Differential/ode_params.jl")
include("../CommonLibrary/struct_manipulation.jl")

# Settings
np = 22 # number of parameters to be analysed

# Define the problem object
global problem = create_problem(
    max_day=100.0,
    treatment=CBD_IL_12_ver7,
    model="takuya",
    param_struct=christian
)

_p, _u0 = struct_split(christian) # Default parameters

# Wrap the problem inside a function that solves it with a given set of params 
wrapped_response = function (params_vector)
    # Assume vector = [p; vl0]
    p = params_vector[1:end-1] 
    # println(size(params_vector))
    # println(size(p))  
    u0 = [0.0084; 9.56; 4.95; params_vector[end]; 0]

    # prob1 = remake(problem;p=p, u0=u0)
    sol = solve(problem;p=p, u0=u0, saveat=0.1)
    area = sum((sol[4,:] + sol[5,:]) .* 0.1)

    return area
    # return (sol[4,end] + sol[5,end])
end 

# Lower and upper bounds of the solution 
vector = [_p; _u0[4];]
width_factor = 0.5
lb = (1-width_factor)*vector
ub = (1+width_factor)*vector

# Labels must be in same order as parameter vector
labels = [
    "tem", "td", "tₗ", "NA", "NA",
    "k₁", "k₂", "k₃", "k₄", "k5", "k6",
    "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8",
    "s1", "s2", "vₗ₀"
    ]
res_sens = gsa(wrapped_response, eFAST(),[[lb[i], ub[i]] for i=1:np], samples=2_000)
# print(size(res_sens.S1))
# scatter(1:np, (res_sens.S1+res_sens.ST)'; ticks=(1:np,labels))
# println("n₁ : S1=$(res_sens.S1[1, 23]) | ST=$(res_sens.ST[1, 23])")

groupedbar([res_sens.S1; res_sens.ST]'; 
    bar_position=:stack,
    xticks=(1:np,labels),
    group=repeat(["First-Order", "Total-Order"], inner=np))

