using DifferentialEquations

include("./ode_core.jl")
include("./ode_params.jl")
include("../treatments_lib.jl")

struct updateParams 
    k6::Float64
    d1::Float64
    s2::Float64
end

"""
Function to perform tumour simulations in restricted parameter space.
Pass the parameters of interest in a struct, the unspecified ones will be
inferred from the average values (from GA)
"""
function restricted_simulation(
    restricted,
    tr::Treatment = CBD_IL_12_ver7,
    base::baseParams = christian,
    )

    # Update the default parameters, to be passed to the DDE model
    for name in fieldnames(typeof(restricted))
        val = :($restricted.$name)
        setproperty!(base, name, eval(val))
    end
    
    # Unpack the updated structure into a vector
    vectorised_struct = Vector{Any}(undef, 22)
    for (i, name) in enumerate(fieldnames(typeof(base)))
        vectorised_struct[i] = eval(:($base.$name))
    end
    p = [vectorised_struct[1:end-1]; tr]
    u0 = [vectorised_struct[end]; 0]
    
    t_span = (0.0, 27.0)
    h(p, t; idxs::Int) = 0.0

    prob_dde = DDEProblem(full_immune_response, u0, h, t_span, p)
    return prob_dde
end