# File with various tools and distributions not offered by standard packages.
using Random
using Distributions

# TODO: vectorize the function, slow for high n
"""
Julia implementation of the binormal distribution.
"""
function binorm(alpha=0.5, µ1=5, µ2=10, s1=1, s2=1,n =1)
    # Output placeholder
    sampled_numbers = Array{Float64}

    # Declartion of the Distribution objects
    norm1 = Normal(µ1, s1)
    norm2 = Normal(µ2, s2)
    
    for (id, outcome) in enumerate(rand(n))
        if outcome <= alpha 
            sampled_numbers[id] = rand(norm1, 1)
        else
            sampled_numbers[id] = rand(norm2, 1)
        end
    end

    return sampled_numbers
end


