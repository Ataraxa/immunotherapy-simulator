include("../Model/ode_model.jl")

function solve_vector(params)
    p = params[1:21]
    u0 = params[22:end]
    t_span = (0.0, 27.0)
    h(p, t; idxs::Int) = 0.0
    problem = DDEProblem(full_immune_response, u0, h, t_span, p)
    
    sol = solve(problem; saveat=0.1)

    return(sol)
end

placebo = Dict(
    "t_in" => 1,
    "t_in12" => 1,
    "t_inCPI" => 1,
    "active_cbd" =>false,
    "active_il12" =>false,
    "active_cpi" => false,
    "active_days" => [7, 10],
    "mean" => [63.2878, 280.4874],
)

cbd_il_12 = Dict(
    "t_in" => [7],
    "t_in12" => 1,
    "t_inCPI" => 1,
    "active_cbd" =>true,
    "active_il12" =>false,
    "active_cpi" => false,
    "active_days" => [7,10,12,14,16,18,20,25,27],
    "mean" => [62.1504, 102.2926, 99.4369, 54.1510, 32.8043, 11.0508, 5.3523, 12.7386, 41.0360]
)

cbd_il_9_14 = Dict(
    "t_in" => [9, 17],
    "t_in12" => 1,
    "t_inCPI" => 1,
    "active_cbd" =>true,
    "active_il12" =>false,
    "active_cpi" => false,
    "active_days" => [8, 11, 13, 15, 17, 19, 21, 23, 25],
    "mean" => [81.6258, 192.2204, 185.5505, 151.3828, 114.7650, 93.7092, 92.8664, 140.6948, 224.7854]
)

il_12_7 = Dict(
    "t_in" => 1,
    "t_in12" => 7,
    "t_inCPI" => 1,
    "active_cbd" =>false,
    "active_il12" =>true,
    "active_cpi" => false,
    "active_days" => [6, 11, 14],
    "mean" => [64.1556, 274.8722, 366.7075],
)

cpi_9_14 = Dict(
    "t_in" => 1,
    "t_in12" => 1,
    "t_inCPI" => 9,
    "active_cbd" =>false,
    "active_il12" =>false,
    "active_cpi" => true,
    "active_days" => [8, 11, 13],
    "mean" => [81.7531, 248.9018, 497.5617],
)

combo_therapy = Dict(
    "t_in" => [9, 14],
    "t_in12" => 1,
    "t_inCPI" => 9,
    "active_cbd" =>true,
    "active_il12" =>false,
    "active_cpi" => true,
    "active_days" => [8,11,13,15,17,19,21,23,25],
    "mean" => [77.3493,207.7220,174.8117,142.1405,102.3850,83.0346,60.0445,47.6723,49.3991]
)

treatments_available = Dict(
    "placebo" => placebo,
    "cbd_7" => cbd_il_12,
    "cbd_914" => cbd_il_9_14,
    "il_7" => il_12_7,
    "cpi" => cpi_9_14,
    "combo" => combo_therapy)