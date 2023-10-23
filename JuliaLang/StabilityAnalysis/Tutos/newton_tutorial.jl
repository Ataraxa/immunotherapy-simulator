using BifurcationKit
using ForwardDiff # for auto-differentiation
using LinearAlgebra
using Plots: plot, scatter

function mixed_jacobian(z, k, f, J)
    x = z[1:end-1]; p = z[end]
    # start creating the mixed space jacobian
    j = J(x, p)
    # to the state space jacobian add one more column, derivative towards p
    pder = ForwardDiff.derivative(p -> f(x, p), p)
    Jmixed = hcat(j, pder)
    # add the last row, which is 1 for the `k` entry, 0 everywhere else
    last_row = zeros(length(z)); last_row[k] = 1.0
    Jfinal = vcat(Jmixed, last_row')
    return Jfinal
end

function newton_step!(zⱼ, zpred, i, f, J, δ)
    Jfinal = mixed_jacobian(zⱼ, i, f, J)
    xⱼ = zⱼ[1:end-1]; pⱼ = zⱼ[end]
    g = f(xⱼ, pⱼ)
    gz = vcat(g, zⱼ[i] - zpred[i])
    zⱼ₊₁ = zⱼ - Jfinal \ gz
    return zⱼ₊₁
end

function corrector(zpred, f, J; δ = 0.9, max_steps = 200, ε = 1e-6, k = 1)
    c = 0
    zⱼ = zpred
    zⱼ₊₁ = newton_step!(zⱼ, zpred, k, f, J, δ)
    while norm(zⱼ₊₁ - zⱼ) > ε
        zⱼ = zⱼ₊₁
        zⱼ₊₁ = newton_step!(zⱼ, zpred, k, f, J, δ)
        c += 1
        if c > max_steps
            @warn("Newton did not converge.")
            return (zⱼ₊₁, false)
        end
    end
    return zⱼ₊₁, true
end

function predictor(zs, dz0)
    if length(zs) == 1
        return zs[end]
    elseif length(zs) == 2 # 1 entry is z0, 2nd entry is 1st found fixed point
        return zs[end] .+ dz0
    else
        return 2zs[end] .- zs[end-1]
    end
end

function continuation!(zs, f, J; dz0, pmin, pmax)
    zpred = predictor(zs, dz0)
    (pmin ≤ zpred[end] ≤ pmax) || return false
    zˣ, success = corrector(zpred, f, J)
    push!(zs, zˣ)
    return success
end

function continuation(f, J, x0, p0; pmin, pmax, dp0, dx0, N = 1000)

    z0 = vcat(x0, p0); zs = [z0]; dz0 = vcat(dx0, dp0)

    ps = [p0]
    xs = [x0]
    stability = Bool[]
    for i in 1:N
        success = continuation!(zs, f, J; dz0, pmin, pmax)
        # Stop iteration if we exceed given parameter margins
        success || break
        # Detect stability of found fixed point (needs `Array` coz of StaticArrays.jl)
        eigenvalues = eigvals(J(zs[end][1:end-1], zs[end][end]))
        isstable = maximum(real, eigenvalues) < 0
        push!(stability, isstable)
    end
    xs = [z[1:end-1] for z in zs]
    ps = [z[end] for z in zs]
    popfirst!(xs); popfirst!(ps) # remove initial guess
    return xs, ps, stability
end

# function maasch_rule(u, p)
#     x, y, z = u
#     q, r, s = 1.2, 0.8, 0.8
#     dx = -x - y
#     dy = -p*z + r*y + s*z^2 - z^2*y
#     dz = -q*(x + z)
#     return [dx, dy, dz]
# end

# function maasch_jacob(u, p)
#     x, y, z = u
#     q, r, s = 1.2, 0.8, 0.8
#     return [-1     -1  0;
#             0   (r - z^2)  (-p + 2z*s - 2z*y);
#             -q   0   -q]
# end

function maasch_rule(u, p)
    x = u
    dx = p .+ x.^2
    return [dx]
end

function maasch_jacob(u, p)
    x = u
    return [2*x]
end


pmin = -1.0
pmax = 1.0
δ = 0.9
p0 = 0.9
x0 = [0]
dp0 = 0.1
dx0 = [0.01]

xs, ps, stability = continuation(maasch_rule, maasch_jacob, x0, p0;
    pmin, pmax, dp0, dx0
)

colors = [s ? :blue : :red for s in stability]
p = scatter(ps, [x[1] for x in xs]; color = colors, markerstrokecolor = colors, xlabel = "p", ylabel = "x", label = "")