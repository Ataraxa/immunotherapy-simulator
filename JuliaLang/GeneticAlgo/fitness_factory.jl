using Match
include("../Model/mechanistic_model.jl")

function create_fitness(; model="takuya")

    initial_prob = create_problem(model=model)
    shift = (model == "w/feedback")

    "Argument format = u0 DOES NOT contain the 0"
    function fitness(params)
        ### Variable initialisation
        prob = initial_prob
        p = params[1:21+shift]
        u0 = params[(22+shift):end]
        debug_data = Vector{Vector{Tuple}}(undef, 6) # this is simply for debug purpose - will remove later
        error_per_treatment = zeros(6)
        j = 1

        ### Iterate over all treatments
        for treatment in treatments_available
            debug_data[j] = Vector{Tuple}(undef, 0)

            sol = solve(prob; p=p, u0=u0, saveat=0.1)
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

    return fitness
end