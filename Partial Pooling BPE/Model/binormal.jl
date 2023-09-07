# TODO - clean up the code
using Distributions

struct BiNormal{T<:Real} <: Distribution{Univariate,Continuous}
    a
    # b::Turing.ForwardDiff.Dual{Turing.ForwardDiff.Tag{Turing.TuringTag, Float64}, Float64} 
    b
    c
    d 
    e::T

            # inner constructor function to instantiate new BiNormal objects
    function BiNormal{T}(a, b, c, d, e::T; check_args = true) where {T<:Real}
            check_args && Distributions.@check_args(BiNormal, c>0 && d>0)
        return new{T}(a,b,c,d,e)
    end
end

# constructor functions for implicitly supplied type
# constructor for no type and params Float64
        function BiNormal(a::Float64,b::Float64,c::Float64,d::Float64,e::Float64; check_args = true)
    return BiNormal{Float64}(a,b,c,d,e,check_args = check_args)
end


# constructor for real params - use promote to make aprams the same type
BiNormal(a::Real, b::Real, c::Real, d::Real, e::Real) = BiNormal(promote(a,b,c,d,e)...)
BiNormal(a::Integer, b::Integer, c,d::Integer,e::Float64) = 
    BiNormal(float(a),float(b), float(c),float(d),float(e))

######################END NEW TYPE SECTION

##### BEGIN EIGHT METHODS

import Base.rand, StatsBase.params
import Random, Distributions, Statistics, StatsBase
using Random
# 0 - helper function
StatsBase.params(d::BiNormal) = (d.a, d.b, d.c, d.d, d.e)

#1 rand(::AbstractRNG, d::UnivariateDistribution)
function Base.rand(rng::AbstractRNG, d::BiNormal)
    (µ1, µ2, s1, s2, a) = StatsBase.params(d)
    norm1 = Normal(µ1, s1)
    norm2 = Normal(µ2, s2)

    return (rand() < a) ? rand(norm1) : rand(norm2)
end


#2 sampler(d::Distribution) - works for sampler(rng::AbstractSampler, d::Distribution)
Distributions.sampler(rng::AbstractRNG,d::BiNormal) = Base.rand(rng::AbstractRNG, d::BiNormal)


#3 logpdf(d::UnivariateDistribution, x::Real)
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

#4 cdf(d::UnivariateDistribution, x::Real)
function Distributions.cdf(d::BiNormal{T}, x::Real) where T<:Real
    (a, b) = params(d)
    if x <= 0
        return(zero(T)) ## equivalent of zero for type T
    elseif x >= 1
        return(one(T)) ## equivalent of 1 for type T
    else
        return(1 - (1-x^a)^b)
    end
end


#5 quantile(d::UnivariateDistribution, q::Real)
function Statistics.quantile(d::BiNormal{T}, x::Real) where T<:Real
    (a, b) = params(d)
    if x <= 0
        return(zero(T)) ## equivalent of zero for type T
    elseif x >= 1
        return(one(T)) ## equivalent of 1 for type T
    else
        return((1 - (1-x)^(1/b))^(1/a))
    end
end


#6 minimum(d::UnivariateDistribution)
function Base.minimum(d::BiNormal)
    return(0)
end


#7 maximum(d::UnivariateDistribution)
function Base.maximum(d::BiNormal)
    return(+Inf)
end


#8 insupport(d::UnivariateDistribution, x::Real)
function Distributions.insupport(d::BiNormal)
    insupport(d::BiNormal, x::Real) = zero(x) <= x <= +Inf
end

######################END Eight Methods
