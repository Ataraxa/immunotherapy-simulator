"""
Set of tests to check that the custom binormal distribution is correct according
to the definition given in the proposal.
"""

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
binormal = BiNormal{Float64}(3.0, 12.0, 1.0, 1.0, 0.5)
norm1 = Normal(5, 1)
check_pdf_bruteforce(binormal)
return 0