using JLD2
using Plots: plot, plot!, scatter!, RGB
using Evolutionary

include("../Model/ode_model.jl")
include("../Library/optimisation_lib.jl")
params = load_object("Res/ga_res.jld2")
what_to_plot = "benchmark"

# Paramer vector
if what_to_plot == "ga_res"
    params = Evolutionary.minimizer(params)
    global opt_p = params[1:21]
    global _, true_p = get_default_values() 
    global u0 = [params[22:end]; 0]
elseif what_to_plot == "benchmark"
    u0, opt_p = get_default_values()
end

# Default value for DDE problem
t_span = (0.0, 27.0)
h(p, t; idxs::Int) = 0.0

# Solve and plot for each treatment
function sotr(i)
    tr = treatments_available[i]
    passed_params = [opt_p; tr]
    prob_dde = DDEProblem(full_immune_response, u0, h, t_span, passed_params)
    sol = solve(prob_dde; saveat=0.1)

    return sol[4,:] + sol[5,:]
end


layout = plot(layout=(2,3))
titles=["Placebo", "CBD-IL-12 (Day 7)", "DBC-IL-12 (Day 9 & 14)", "IL-12 (Day 7)", "CPI (Day 9 & 14)", "CPI + CBD-IL-12"]
for i = 1:6
    simulated = sotr(i)
    tr = treatments_available[i]
    xaxis = tr["active_days"]
    plot!(layout, 0.0:0.1:27.0,  simulated;
    
    linecolor=RGB(0, 0.4470, 0.7410),
    titlefontsize=7,
    xguidefontsize=5,
    yguidefontsize=5,
    legend=false, 
    subplot=i, 
    title=titles[i], 
    xlabel="Time (days)",
    ylabel="Tumour volume (mm³)")

    scatter!(layout, xaxis, simulated[xaxis*10];
        subplot=i, mc=RGB(0, 0.4470, 0.7410), ms=:2)

    plot!(layout, xaxis, tr["mean"];
        subplot=i, linecolor=RGB(0.8500, 0.3250, 0.0980))
        
    scatter!(
        layout, xaxis, tr["mean"];
        subplot=i, mc=RGB(0.8500, 0.3250, 0.0980), ms=2)
end


plot(layout)
