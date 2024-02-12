using JLD2
using Plots
include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")

abc_results = load_object("Results/abc/ABC-local-0.jld2")
popu = abc_results.population[7]
prob = create_problem(model="takuya")
params = copy(christian_true_params)

layout = plot(layout=(2,2))

for row = 1:size(popu)[1]
    params[[11,12,21]] .= popu[row,:] 
    pred = solve(prob; p=params, saveat=0.1)
    v = pred[4,:] + pred[5,:]

    combined_pred = vcat(pred[1:3,:], reshape(v, 1, length(v)))
    series = ["IFNγ", "CD8+", "PD-1", "Tumour volume"]

    for i in 1:4
        plot!(layout, pred.t, combined_pred[i,:];
        sp=i,
        xlabel="Time (day)",
        title=series[i],
        frame=:box,
        legend=false)
        if i == 4
            ylims!((0, 600); sp=i)
        end
    end

end
 display(layout)