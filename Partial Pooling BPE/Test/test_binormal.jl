"""
Set of tests to check that the custom binormal distribution is correct according
to the definition given in the proposal (see Documents).
"""

using Turing
using Distributions
using StatsPlots: plot, scatter
using StatsBase: sample
using HDF5

include("../Model/binormal.jl")

"Plots the pdf of target distribution by naive sampling"
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

"Plots pdf of target distribution using MCMC sampling"
function coin_flip_check(distrib::Distribution)
    @model function coinFlip()
        θ ~ distrib 
    end

    chain = sample(coinFlip(), NUTS(0.7), 100)
    return(chain)
end

"Plots pdf of hyper-binormal distribution"
function hyper_coin_check()
    @model function hyperCoinFip()
        µ1 ~ Normal{Float64}(7.0, 1.0)
        µ2 ~ Normal{Float64}(2.0, 1.0)
        a ~ Beta(1,1)
        s ~ truncated(Normal(0, 1); lower=1)
        θ ~ BiNormal{Float64}(µ2, µ1, s, s, a)
    end

    chain = sample(hyperCoinFip(), NUTS(0.7), MCMCSerial(), 1000, 3)
    return(chain)
end

# Define distributions to test or act as positiv control
binormal = BiNormal{Float64}(0.2, 0.8, 0.1, 0.1, 0.5)
normal = Normal(5, 1)
beta_uniform = Beta(2, 2)

# Perform check 
res = coin_flip_check(normal)
# res = hyper_coin_check()
# plot(res)

h5open("test.h5", "w") do f 
    write(f, res)
end