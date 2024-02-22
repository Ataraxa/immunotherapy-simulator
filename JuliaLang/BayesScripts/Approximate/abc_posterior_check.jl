using JLD2
using Plots
include("../../Model/mechanistic_model.jl")
include("../../Model/Bayesian/priors.jl")
plotlyjs()
abc_results = load_object("Results/abc/ABC-local-2.jld2")
popu = abc_results.population[end]
prob = create_problem(model="takuya")
params = copy(christian_true_params)

# varidx = collect(11:21)
varidx = [11, 12, 21]
layout = plot(layout=(2,2))
N = size(popu)[1]
for row = 1:N
    params[varidx] .= exp.(popu[row,:])
    # println.(popu[row,:])
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
        legend=false,
        color=:grey,
        alpha=0.5
        )
        # if i == 4
        #     ylims!((0, 600); sp=i)
        # end
    end

end

params[[11,12,21]] .= [1.98814060169, 2.10643371507, 1.65262881505]
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
        ylims!((0, 700); sp=i)
    end
end
display(layout)