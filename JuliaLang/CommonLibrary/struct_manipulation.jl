using UnPack

include("../Model/Differential/ode_params.jl") # For type hinting

## This structure represents the restricted parameter space
mutable struct updateParams1
    k6::Float64
end

struct updateParams3
    k6::Float64
    d1::Float64
    s2::Float64
end

"Function to return (p, u0) from parameter struct"
function struct_split(input_struct::baseParams)
    @unpack t_d, t_delay, t_last, t_delay12, t_last12, 
    k1, k2, k3, k4, k5, k6, 
    d1, d2, d3, d4, d5, d6, d7, d8, 
    s1, s2, u0 = input_struct

    p = [t_d, t_delay, t_last, t_delay12, t_last12, 
    k1, k2, k3, k4, k5, k6, 
    d1, d2, d3, d4, d5, d6, d7, d8, 
    s1, s2]

    return p, u0
end

"Function to return (p, u0) from parameter struct"
function struct_split(input_struct::neoParams)
    @unpack t_d, t_delay, t_last, t_delay12, t_last12, 
    k1, k2, k3, k4, k5, k6, 
    d1, d2, d3, d4, d5, d6, d7, d8, 
    s1, s2, n1, u0 = input_struct

    p = [t_d, t_delay, t_last, t_delay12, t_last12, 
    k1, k2, k3, k4, k5, k6, 
    d1, d2, d3, d4, d5, d6, d7, d8, 
    s1, s2, n1]

    return p, u0
end