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
include("../Model/template_distrib.jl")
include("../Model/expand_from_template.jl")

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

    chain = sample(coinFlip(), HMC(0.1, 5), 1000)
    return(chain)
end

# Define distributions to test or act as positiv control
binormal = BiNormal{Float64}(3.0, 12.0, 1.0, 1.0, 0.5)
norm1 = Normal(5, 1)
bet23 = Beta(2, 2)
dummy = Kuma3{Int64}(10, 2, 3, 1, 1)

# Perform check 
res = coin_flip_check(dummy)
plot(res)