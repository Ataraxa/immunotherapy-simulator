
using Distributions

struct Kuma2{T<:Real} <: Distribution{Univariate,Continuous}
    a::T  ## first parameter of Kumaraswamy Distribution
    b::T  ## second parameter of Kumraswamy Distribution

    # inner constructor function to instantiate new Kuma2 objects
    function Kuma2{T}(a::T, b::T; check_args = true) where {T<:Real}
        check_args && Distributions.@check_args(Kuma2, a>0 && b>0)
        return new{T}(a,b)
    end
end

# test instantiation
Kuma2{Float64}(2.0,12.0)

# constructor functions for implicitly supplied type
# constructor for no type and params Float64
function Kuma2(a::Float64, b::Float64; check_args = true)
    return Kuma2{Float64}(a,b,check_args = check_args)
end

# constructor for real params - use promote to make aprams the same type
Kuma2(a::Real, b::Real) = Kuma2(promote(a,b)...)
Kuma2(a::Integer, b::Integer) = Kuma2(float(a),float(b))

# test instantiation
Kuma2(2,12)

######################END NEW TYPE SECTION

##### BEGIN EIGHT METHODS

import Base.rand, StatsBase.params
import Random, Distributions, Statistics, StatsBase
using Random
# 0 - helper function
StatsBase.params(d::Kuma2) = (d.a, d.b)

#1 rand(::AbstractRNG, d::UnivariateDistribution)
function Base.rand(rng::AbstractRNG, d::Kuma2)
    (a, b) = params(d)
    u = rand(rng) # get uniform rv for inverse transorm method
    return((1 - (1-u)^(1/b))^(1/a))
end


#2 sampler(d::Distribution) - works for sampler(rng::AbstractSampler, d::Distribution)
Distributions.sampler(rng::AbstractRNG,d::Kuma2) = Base.rand(rng::AbstractRNG, d::Kuma2)


#3 logpdf(d::UnivariateDistribution, x::Real)
function Distributions.pdf(d::Kuma2{T}, x::Real) where {T<:Real}
    (a, b) = params(d)
    if x<= 0 
        return zero(T) ## equivalent of zero for type T
    elseif x>=1
        return zero(T)
    else
        return(a*b*x^(a-1)*(1-x^a)^(b-1))
    end
end

Distributions.logpdf(d::Kuma2, x::Real) = log(pdf(d,x))

#4 cdf(d::UnivariateDistribution, x::Real)
function Distributions.cdf(d::Kuma2{T}, x::Real) where T<:Real
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
function Statistics.quantile(d::Kuma2{T}, x::Real) where T<:Real
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
function Base.minimum(d::Kuma2)
    return(0)
end


#7 maximum(d::UnivariateDistribution)
function Base.maximum(d::Kuma2)
    return(1)
end


#8 insupport(d::UnivariateDistribution, x::Real)
function Distributions.insupport(d::Kuma2)
    insupport(d::Kuma2, x::Real) = zero(x) <= x <= one(x)
end

######################END Eight Methods

#### Test to see if it works

y = Kuma2(2,8)

rand(y,10)


using Turing

# Define a simple coin flip model with unknown prob θ
@model function coinFlip(x)
  θ ~ Kuma2(2,2)  
  # θ ~ Beta(2,2)  
  x ~ Bernoulli(θ)
end

#  Flip 1 coin (heads) - Run sampler, collect results
chn = sample(coinFlip(1), HMC(0.1, 5), 1000)