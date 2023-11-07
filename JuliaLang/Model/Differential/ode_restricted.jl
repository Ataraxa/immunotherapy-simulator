using DifferentialEquations

include("./ode_params.jl")
include("../../CommonLibrary/struct_manipulation.jl")

## This structure represents the restricted parameter space
struct updateParams 
    k6::Float64
    # d1::Float64
    # s2::Float64
end

"""
Use a default parameter structure to complete a restricted-space parameter 
vector

Returns a Tuple (p, u0)
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

    # Return base structure converted to vector
    return struct_split(base)
end

