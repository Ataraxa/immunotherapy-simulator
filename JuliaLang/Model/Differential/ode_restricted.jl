using DifferentialEquations

include("./ode_params.jl")
include("../../CommonLibrary/struct_manipulation.jl")

"""
Use a default parameter structure to complete a restricted-space parameter 
structure.

Returns a Tuple (p, u0)
"""
function repack_params(
    restricted::Union{updateParams1, updateParams3},
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

