# File with various tools and distributions not offered by standard packages.
import Base.rand
import StatsBase
using Random
using Distributions: sampler, logpdf

"""
Binormal distribution, as defined in the analysis proposal, implemented
in Julia.
As for now, it only supports Float64 arguments, but this will be fixed later.
It also enforces truncation at 0, but more flexibility will be added later on.
"""
struct BiNormal{T<:Real} <: Distribution{Univariate, Continuous}
    µ1::T 
    µ2::T 
    s1::T
    s2::T
    a::T
    
    function (BiNormal{T}(µ1::T, µ2::T, s1::T, s2::T, a::T, ; check_args=false) 
            where {T<:Real})
        check_args && Distributions.@check_args(BiNormal, s1>0 && s2>0 && 0<a<1)
        return new{T}(µ1, µ2, s1, s2, a)
    end

end

# Helper function to get the parameters of a distribution
StatsBase.params(d::BiNormal) = (d.µ1, d.µ2, d.s1, d.s2, d.a)

# Reimplementation of a sampling method
function Base.rand(rng::AbstractRNG, d::BiNormal)
    (µ1, µ2, s1, s2, a) = StatsBase.params(d)
    norm1 = Normal(µ1, s1)
    norm2 = Normal(µ2, s2)

    return (rand() < a) ? rand(norm1) : rand(norm2)
end
Distributions.sampler(rng::AbstractRNG,d::BiNormal) = Base.rand(rng::AbstractRNG, d::Kuma)

# Reimplementation of a logpdf method
function Distributions.pdf(d::BiNormal{T}, x::Real) where {T<:Real}
    (µ1, µ2, s1, s2, a) = StatsBase.params(d)
    norm1 = Normal(µ1, s1)
    norm2 = Normal(µ2, s2)

    if x<= 0
        return(zero(T))
    else
        return(a*pdf(norm1,x) + (1-a)*pdf(norm2,x))
    end
end
Distributions.logpdf(d::BiNormal, x::Real) = log(pdf(d,x))

# Main

return "> no error" 


