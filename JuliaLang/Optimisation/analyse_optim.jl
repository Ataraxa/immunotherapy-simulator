using JLD2
using Plots: plot, plot!, scatter!, RGB
using Evolutionary

include("../Model/Differential/ode_model.jl")
include("./opt_lib.jl")

what_to_plot = "benchmark"
filename = "Res/ga_res.jld2"

# Extract paramater vector from files
if what_to_plot == "ga_res"
    # params = load_object(filename)
    params = [1.8728241605057228, 0.4879233369602657, 4.8934754961589695, 
    3.697285769720824, 1.0754112854143516, 0.22193617235217952, 
    6.081963924603686, 74.60149998881774, 929.5659753695654, 5.815052050675676, 
    0.5394946971270578, 10.177238034144443, 10.610875554141064, 
    1.3262458309759877, 4.486099857321717, 0.01602458080336531, 
    0.03416214467083285, 59.611074450040185, 0.5697851077620799, 
    14.068450976326906, 0.4108126097661595, 0.007911968983771049, 
    8.851474387488993, 5.999229998549251, 5.573047884761726]

    # params = Evolutionary.minimizer(params)

    global opt_p = params[1:21]
    global _, true_p = get_default_values() 

    global u0 = [params[22:end]; 0]
elseif what_to_plot == "benchmark"
    u0, opt_p = get_christian()
    # opt_p = (1 - 0.15) .* opt_p
end

# AUX ------------------------------------------------------
# Default value for DDE problem
t_span = (0.0, 27.0)
h(p, t; idxs::Int) = 0.0

# Solve and plot for a given treatment (passed as index)
function sotr(i)
    tr = treatments_available[i]
    passed_params = [opt_p; tr]
    prob_dde = DDEProblem(full_immune_response, u0, h, t_span, passed_params)
    sol = solve(prob_dde; saveat=0.1)

    return sol[4,:] + sol[5,:]
end
# AUX (end) -----------------------------------------------

# Plot tumour evolution or each treatment
layout = plot(layout=(2,3))
titles=["CBD-IL-12 (Day 7)", "CBD-IL-12 (Day 9 & 14)", "CPI + CBD-IL-12", "CPI (Day 9 & 14)", "IL-12 (Day 7)","Placebo"]
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
display(plot(layout))

# Compute fitness 
fitness([opt_p; u0[1:end-1]])





