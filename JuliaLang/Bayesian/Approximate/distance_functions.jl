using UnPack 

include("../../Model/Differential/ode_restricted.jl")
include("../../CommonLibrary/struct_manipulation.jl")

function mse_distance(params, constants, target_data)
    problem, selected_days, step = constants

    params = exp.(params)
    repacked_p, _ = repack_params(updateParams3(params...))
    sim = solve(problem; p=repacked_p, saveat=step)
    sim_data = (sim[4,:] + sim[5,:])[selected_days*trunc(Int, 1/step) .+ 1]

    println(size(sim_data))
    println(size(target_data))
    distance = mean((sim_data - target_data).^2)
    return distance, 1
end