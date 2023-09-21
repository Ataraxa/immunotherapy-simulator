using Evolutionary
using DelimitedFiles

include("../Model/ode_model.jl")

# Fetch known fit parameters for boxConstraints
u0, p = get_default_values()
u0 = u0[1:end-1] # last element is starting dead tumour vol., so irrelevant
target = [p; u0] # 1x25 vector
width_factor = 0.5
lower = [(1 - width_factor) .* target; 0.0] # Add unknown params at the end
upper = [(1 + width_factor) .* target; 1.0]

# Perform optimisation
params = Evolutionary.optimize(
    fitness(),
    BoxConstraints(lower, upper)
    ones(26),
    GA(),
)

function fitness(params)
     for treatment in treatments_available
        sol = solve(params)
        tumour_vol = sol[4,:] + sol[5,:]
        
        i = 1
        error_per_day = zeros()
        for day = treatment["valid_days"]
            _, sol_index = findmin(abs(day .* ones(1, numel(sol.x)) - sol.x))
            in_silico = tumour_vol[sol_index]
            in_vivo = treatment["mean"][i]
            error_per_day[i] = ((in_silico - in_vivo)^2)/1000;
            i += 1
        end

        error_per_treatment(j) = sum(error_per_day)
    end

    error = 0 
    for i = 1:length(error_per_treatment)
        error += error_per_treatment[i]
    end

    return(error)
end

# Save result to file
writedlm("Res/ga_results.csv", params, ',')