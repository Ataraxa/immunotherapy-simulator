"""
Set of tests to check that the custom binormal distribution is correct according
to the definition given in the proposal.
"""

using Turing
using Distributions
using StatsPlots: plot, scatter
using StatsBase: sample

# Custom libraries
include("../Model/binormal.jl")
"""
Plots the pdf of target distribution simply by sampling random numbers from 
it
"""
function check_pdf_bruteforce(d::Distribution; 
        n::Int64=1_000, lower::Int64=0, upper::Int64=15)
    indexed_pdf = zeros(n)
    sampled_pdf = zeros(n)
    sampler = Uniform(lower, upper)

    for i = 1:n
        rndnum = rand(sampler)
        sampled_pdf[i] = pdf(d, rndnum)
        indexed_pdf[i] = rndnum
    end
    display(scatter(indexed_pdf, sampled_pdf))
end

# Main 

function coin_flip_check(distrib::Distribution)
    @model function coinFlip()
        θ ~ distrib 
    end

    chain = sample(coinFlip(), NUTS(0.7), MCMCSerial(), 1000, 3)
    return(chain)
end

function hyper_coin_check()
    @model function hyperCoinFip()
        µ ~ Normal{Float64}(10.0, 1.0)
        θ ~ BiNormal{Float64}(5.0, µ, 0.1, 0.1, 0.5)
    end

    chain = sample(hyperCoinFip(), NUTS(0.7), MCMCSerial(), 1000, 3)
    return(chain)
end

# Define distributions to test or act as positiv control
# binormal = BiNormal{Float64}(0.2, 0.8, 0.1, 0.1, 0.5)
normal = Normal(5, 1)
beta_uniform = Beta(2, 2)

# Perform check 
# res = coin_flip_check(binormal)
res = hyper_coin_check()
plot(res)