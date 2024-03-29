# File with various tools and distributions not offered by standard packages.
import Base.rand
import StatsBase
using Random
using Distributions
using Distributions: sampler, logpdf

"""
Binormal distribution, as defined in the analysis proposal, implemented
in Julia.
As for now, it only supports Float64 arguments, but this will be fixed later.
It also enforces truncation at 0, but more flexibility will be added later on.
"""
struct BiNormal{T<:Real} <: ContinuousUnivariateDistribution
    µ1::T 
    µ2::T 
    s1::T
    s2::T
    a::T
    BiNormal{T}(µ1::T, µ2::T, s1::T, s2::T, a::T) where {T} = 
        new{T}(µ1, µ2, s1, s2, a)
end

function BiNormal(µ1::T, µ2::T, s1::T, s2::T, a::T; check_args=true) where {T<: Real}
    Distributions.@check_args BiNormal (µ1, µ1 > 0)
    return BiNormal{T}(µ1, µ2, s1, s2, a)
end

# Distributions.@distr_support BiNormal 5.0 10.0 1.0 1.0 0.5

### Conversions

### Parameters
StatsBase.params(d::BiNormal) = (d.µ1, d.µ2, d.s1, d.s2, d.a)

### Sampling 

function Base.rand(rng::AbstractRNG, d::BiNormal)
    (µ1, µ2, s1, s2, a) = StatsBase.params(d)
    norm1 = Normal(µ1, s1)
    norm2 = Normal(µ2, s2)

    return (rand() < a) ? rand(norm1) : rand(norm2)
end
Distributions.sampler(rng::AbstractRNG,d::BiNormal) = Base.rand(rng::AbstractRNG, d::Kuma)
Distributions.sampler(d::BiNormal) = Base.rand(d::Kuma)

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


