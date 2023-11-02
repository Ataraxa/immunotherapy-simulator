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
function repack_params(
    restricted::updateParams,
    base::baseParams = christian,
    )

    # Update the default parameters, to be passed to the DDE model
    for name in fieldnames(typeof(restricted))
        val = :($restricted.$name)
        setproperty!(base, name, eval(val))
    end
    
    return prob_dde
end

