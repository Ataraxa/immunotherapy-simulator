using ApproxBayes
using Distributions
using Random
function normaldist(params, constants, targetdata)
    simdata = rand(Normal(params...), 1000)
    ApproxBayes.ksdist(simdata, targetdata), 1
end

Random.seed!(1)
p1 = 2.0
p2 = 0.4
targetdata = rand(Normal(p1, p2), 1000)

setup = ABCRejection(normaldist, #simulation function
  2, # number of parameters
  0.1, #target ϵ
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]); # Prior for each of the parameters
  maxiterations = 10^6, #Maximum number of iterations before the algorithm terminates
  )

# run ABC inference
rejection = runabc(setup, targetdata)