include("./ode_core.jl")
include("./ode_params.jl")
include("../treatments_lib.jl")

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
    for name in fieldnames(typeof(restricted))
        println(typeof(restricted))
        val = getFieldValue(restricted, name)
        setproperty!(base, field, val)
    end

    vectorised_struct = Vector{Float64}(undef, 21)
    for (i, name) in enumerate(fieldnames(typeof(base)))
        p[i] = getFieldValue(base, name)
    end
    p = [vectorised_struct[1:end-1]; tr]
    u0 = vectorised_struct[end]
    
    t_span = (0.0, 27.0)
    h(p, t; idxs::Int) = 0.0

    prob_dde = DDEProblem(full_immune_response, u0, h, t_span, p)
    return prob_dde
end

function getFieldValue(paramStruct, name::Symbol)
    field = Symbol(name)
    code = quote
        (obj) -> obj.$field
    end
    val = eval(code)(paramStruct)
    return val
end