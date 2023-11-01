using DifferentialEquations
using Plots: plot, plot!
include("../Model/Differential/ode_core.jl")

function solve_for_treatment(params, treatment_spec)
    p = [params[1:21]; treatment_spec]
    u0 = [params[22:end]; 0] # Add back the 0 here!

    t_span = (0.0, 27.0)
    h(p, t; idxs::Int) = 0.0

    # println("")
    # println("_______________________________________________")
    # println(u0)

    problem = DDEProblem(full_immune_response, u0, h, t_span, p)
    
    sol = solve(problem; saveat=0.01)

    return(sol)
end

"Argument format = u0 DOES NOT contain the 0"
# function fitness(params)
#     error_per_treatment = zeros(6)
#     j = 1
#     for treatment in treatments_available
#         params = [params; 0]
#         sol = solve_for_treatment(params, treatment)
#         tumour_vol = sol[4,:] + sol[5,:]

#         i = 1
#         error_per_day = zeros(length(treatment["active_days"]))
#         for day = treatment["active_days"]
#             _, sol_index = findmin(broadcast(abs, day .* ones(length(sol.t)) - sol.t))
#             in_silico = tumour_vol[sol_index]
#             in_vivo = treatment["mean"][i]
#             error_per_day[i] = ((in_silico - in_vivo)^2)/1000;
#             i += 1
#         end

#         error_per_treatment[j] = sum(error_per_day)
#         j+=1
#     end

#     error = 0 
#     for i in eachindex(error_per_treatment)
#         error += error_per_treatment[i]
#     end

#     return(error)
# end
function fitness(params, is_debug=false)
    debug_data = Vector{Vector{Tuple}}(undef, 6) # this is simply for debug purpose - will remove later
    error_per_treatment = zeros(6)
    j = 1

    # Iterate over all treatments
    for treatment in treatments_available
        debug_data[j] = Vector{Tuple}(undef, 0)

        sol = solve_for_treatment(params, treatment)
        tumour_vol = sol[4,:] + sol[5,:]
        # plot(treatment["active_days"], treatment["mean"])
        # display(plot!(0.0:0.01:27.0, tumour_vol))
        

        i = 1
        error_per_day = zeros(length(treatment.active_days))
        for day = treatment.active_days
            _, sol_index = findmin(broadcast(abs, day .* ones(length(sol.t)) - sol.t))
            in_silico = tumour_vol[sol_index]
            in_vivo = treatment.mean[i]
            error_per_day[i] = ((in_silico - in_vivo)^2)/1000;

            push!(debug_data[j], (sol_index, in_silico, in_vivo))
            i += 1
        end

        error_per_treatment[j] = sum(error_per_day)
        j+=1
    end

    error = 0 
    for i in eachindex(error_per_treatment)
        error += error_per_treatment[i]*length(error_per_treatment)
    end
    
    println("_______________________________________________")
    println(error)
    println(params)
    
    return(error)
end

# placebo = Dict(
#     "t_in" => 1,
#     "t_in12" => 1,
#     "t_inCPI" => 1,
#     "active_cbd" =>false,
#     "active_il12" =>false,
#     "active_cpi" => false,
#     "active_days" => [7, 10],
#     "mean" => [63.2878, 280.4874],
# )

# cbd_il_12 = Dict(
#     "t_in" => [7],
#     "t_in12" => 1,
#     "t_inCPI" => 1,
#     "active_cbd" =>true,
#     "active_il12" =>false,
#     "active_cpi" => false,
#     "active_days" => [7,10,12,14,16,18,20,25,27],
#     "mean" => [62.1504, 102.2926, 99.4369, 54.1510, 32.8043, 11.0508, 5.3523, 12.7386, 41.0360]
# )

# cbd_il_9_14 = Dict(
#     "t_in" => [9, 14],
#     "t_in12" => 1,
#     "t_inCPI" => 1,
#     "active_cbd" =>true,
#     "active_il12" =>false,
#     "active_cpi" => false,
#     "active_days" => [8, 11, 13, 15, 17, 19, 21, 23, 25],
#     "mean" => [81.6258, 192.2204, 185.5505, 151.3828, 114.7650, 93.7092, 92.8664, 140.6948, 224.7854]
# )

# il_12_7 = Dict(
#     "t_in" => 1,
#     "t_in12" => 7,
#     "t_inCPI" => 1,
#     "active_cbd" =>false,
#     "active_il12" =>true,
#     "active_cpi" => false,
#     "active_days" => [6, 11, 14],
#     "mean" => [64.1556, 274.8722, 366.7075],
# )

# cpi_9_14 = Dict(
#     "t_in" => 1,
#     "t_in12" => 1,
#     "t_inCPI" => 9,
#     "active_cbd" =>false,
#     "active_il12" =>false,
#     "active_cpi" => true,
#     "active_days" => [8, 11, 13],
#     "mean" => [81.7531, 248.9018, 497.5617],
# )

# combo_therapy = Dict(
#     "t_in" => [9, 14],
#     "t_in12" => 1,
#     "t_inCPI" => 9,
#     "active_cbd" =>true,
#     "active_il12" =>false,
#     "active_cpi" => true,
#     "active_days" => [8,11,13,15,17,19,21,23,25],
#     "mean" => [77.3493,207.7220,174.8117,142.1405,102.3850,83.0346,60.0445,47.6723,49.3991]
# )

# treatments_available = Dict(
#     "placebo" => placebo,
#     "cbd_7" => cbd_il_12,
#     "cbd_914" => cbd_il_9_14,
#     "il_7" => il_12_7,
#     "cpi" => cpi_9_14,
#     "combo" => combo_therapy)

# treatments_available = [
#     cbd_il_12, cbd_il_9_14, combo_therapy,
#     cpi_9_14 ,  il_12_7   , placebo
#     ]
