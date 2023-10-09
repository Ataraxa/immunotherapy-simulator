using Plots: plot, plot!,density
using HDF5
using MCMCChains
using MCMCChainsStorage
using Distributions
using StatsPlots
using LaTeXStrings

function draw_population_posterior(chain, n_samples=10_000)
    density_vector_k6 = Vector{Float64}(undef, n_samples)
    density_vector_d1 = Vector{Float64}(undef, n_samples)
    density_vector_s2 = Vector{Float64}(undef, n_samples)
    
    # Sample µ and σ
    post_vec = sample(chain[[:µ_k6, :µ_d1, :µ_s2, :σ_k6, :σ_d1, :σ_s2]], 
        n_samples; replace=true)

    # Sample from population N(µ,σ)
    i = 1
    for p in eachrow(Array(post_vec)) # ie for each sample 
        density_vector_d1[i] = rand(Normal(p[1], p[4]))
        density_vector_k6[i] = rand(Normal(p[2], p[5]))
        density_vector_s2[i] = rand(Normal(p[3], p[6]))
        i += 1
    end
    # plot(density(density_vector_d1); xlim=(0,40))
    # println(median(density_vector_d1))

    plot_font = "Computer Modern"
    default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)

    plt = plot(
        density(filter!(e -> 0 < e < 40, density_vector_d1));
        label="fit",
        xlabel=L"d_1",
        xlim=(0,10))
    display(plot!(truncated(Normal(log(11), 1); lower=0, upper=10); label="true"))
    savefig(plt, "Misc/Images/d1_popu.png")

    plt = plot(
        density(density_vector_k6);
        label="fit",
        xlabel=L"k_6",
        xlim=(-4,2.5))
    display(plot!(Normal(log(0.5), 1); label="true"))
    savefig(plt, "Misc/Images/k6_popu.png")

    plt = plot(
        density(density_vector_s2);
        label="fit",
        xlabel=L"s_2",
        xlim=(-4,2.5))
    plot!(Normal(log(0.3), 1); label="true")
    savefig(plt, "Misc/Images/s2_popu.png")

end

if abspath(PROGRAM_FILE) != @__FILE__
    chain = h5open("Res/hpc-validation_chain-8.h5", "r") do f
        read(f, Chains)
    end

    a=draw_population_posterior(chain)
    # savefig(a, "d1_popu.eps")
end