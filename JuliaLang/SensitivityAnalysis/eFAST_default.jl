using GlobalSensitivity
using DifferentialEquations
using Distributions
using Plots
using StatsPlots: groupedbar
using Match
using Plots.PlotMeasures
include("../Model/treatments_lib.jl")
include("../Model/mechanistic_model.jl")
include("../Model/Bayesian/priors.jl")
# Settings
np = 25 # number of parameters to be analysed

# Define the problem object
global problem = create_problem(
    max_day=100.0,
    treatment=CBD_IL_12_ver7,
    model="takuya",
    params=christian_true_params
)

# Function wrapper that outputs scalar metric from numerical approx
function scalar_metric(num_approx, time_step, metric_type)
    # # Assume vector = [p; vl0]
    # p = params_vector[1:end-1] 
    # # println(size(params_vector))
    # # println(size(p))  
    # u0 = [0.0084; 9.56; 4.95; params_vector[end]; 0]

    # # prob1 = remake(problem;p=p, u0=u0)
    # sol = solve(problem;p=p, u0=u0, saveat=0.1)
    # area = sum((sol[4,:] + sol[5,:]) .* 0.1)
    @match metric_type begin 
        "area" => global output = sum(num_approx .* time_step) 
    end

    return output
end 

function wrapper(param_vector)
    new_p = param_vector[1:21] 
    new_u0 = [param_vector[22:25]; 0]

    sim = solve(problem; p=new_p, u0=new_u0, saveat=0.1)
    tumour = sim[4,:] + sim[5,:]

    return sum(tumour .* 0.1)
end

# Lower and upper bounds of the solution 
vector = christian_true_params
width_factor = 0.5
lb = (1-width_factor)*vector
ub = (1+width_factor)*vector

# Labels must be in same order as parameter vector
labels = [
    "t_d", "t_delay", "t_last", "t_delay12", "t_last12",
    "k1", "k2", "k3", "k4", "k5", "k6",
    "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8",
    "s1", "s2", "g0", "c0", "p0", "v0"
    ]
res_sens = gsa(wrapper, eFAST(),[[lb[i], ub[i]] for i=1:np], samples=2000)
# print(size(res_sens.S1))
# scatter(1:np, (res_sens.S1+res_sens.ST)'; ticks=(1:np,labels))
# println("n₁ : S1=$(res_sens.S1[1, 23]) | ST=$(res_sens.ST[1, 23])")

groupedbar([res_sens.ST;res_sens.S1]'; 
    bar_position=:stack,
    xrotation=90,
    xticks=(1:np,labels),
    xlabel="Model Parameters",
    bottom_margin=7mm,
    dpi=1000,
    group=repeat(["Total-Order","First-Order" ], inner=np)
)

savefig("proutiprouta.png")



