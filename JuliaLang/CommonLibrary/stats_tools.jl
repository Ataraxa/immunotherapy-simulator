using Distributions

function confidence_interval(vector, level)
    µ = mean(vector)
    σ = std(vector)
    n = length(vector)

    d = TDist(n-1)
    mult = quantile(d, (1-level)/2)
    lb = µ - mult * σ/sqrt(n)
    ub = µ + mult * σ/sqrt(n)
    return [lb, up]
end