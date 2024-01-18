using Distributions
using ApproxBayes
using Random
Random.seed!(1)

# Simulator (params in Vect{Float64})
function normaldist(params, constants, targetdata)
    println(typeof(params))
    simdata = rand(Normal(params...), 1000)
    ApproxBayes.ksdist(simdata, targetdata), 1 # need to return a Tuple
end

# Target data
p1 = 2.0
p2 = 0.4
targetdata = rand(Normal(p1, p2), 1000)

# Inference
setup = ABCSMC(normaldist, #simulation function
  2, # number of parameters
  0.1, #target ϵ
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]), #Prior for each parameter
  )

smc = runabc(setup, targetdata, verbose = true, progress = true)