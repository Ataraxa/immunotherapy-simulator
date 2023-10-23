using Plots
using BifurcationKit

Fbp(u, p) = @. u * (1 - p.μ - u / 600)


# bifurcation problem
prob = BifurcationProblem(Fbp, [200.], (μ = 0.,),
	# specify the continuation parameter
	(@lens _.μ),
	record_from_solution = (x, p) -> x[1])

# options for continuation
opts_br = ContinuationPar(
	# parameter interval
	p_max = 2.0, p_min = 0.0,
	# detect bifurcations with bisection method
	# we increase the precision of the bisection
	n_inversion = 4)

# automatic bifurcation diagram computation
diagram = bifurcationdiagram(prob, PALC(),
	# very important parameter. This specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches
	# at most in the present case.
	2,
	(args...) -> opts_br,
	)

plot(diagram)