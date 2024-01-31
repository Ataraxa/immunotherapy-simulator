z = randn(10000)

@model dottilde(y) = begin
    σ = Vector{Float64}(undef, 2)
    μ ~ Normal(0,1)
    σ::Vector{Float64} .~ Gamma(2,1)
    y .~ Normal(μ,σ)
end

q = vi(dottilde(z), ADVI(1,2000))