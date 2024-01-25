using Plots: plot, plot!,density, vline!
using HDF5
using MCMCChains
using MCMCChainsStorage
using Distributions
using StatsPlots
using LaTeXStrings

function draw_population_posterior(chain, n_samples=1000)
    density_vector_k6 = Vector{Float64}(undef, n_samples)
    density_vector_d1 = Vector{Float64}(undef, n_samples)
    density_vector_s2 = Vector{Float64}(undef, n_samples)

    k6 = namesingroup(chain, :k6)
    d1 = namesingroup(chain, :d1)
    s2 = namesingroup(chain, :s2)
    
    # Sample µ and σ
    i = 1
    post_vec = sample(chain[[k6[1], d1[1], s2[1]]], n_samples; replace=false)
    for p in eachrow(Array(post_vec))
        density_vector_k6[i] = p[2]
        density_vector_d1[i] = p[1]
        density_vector_s2[i] = p[3]
        i+=1
    end
    # a = vline([log(11)], linewidth=10)
    plt = plot(density(density_vector_d1); legend=false)
    display(plot!([log(11)], seriestype="vline"; xlabel=L"d_1"))
    savefig(plt, "Misc/Images/d1_1.pdf")

    plt = plot(density(density_vector_k6); legend=false)
    display(plot!([log(0.46)], seriestype="vline"; xlabel=L"k_6"))
    savefig(plt, "Misc/Images/k6_1.pdf")

    plt = plot(density(density_vector_s2); legend=false)
    display(plot!([log(0.16)], seriestype="vline"; xlabel=L"s_2"))
    savefig(plt, "Misc/Images/s2_1.pdf")

    # display(plot(density(density_vector_k6); title=L"k_6"))
    # display(plot(density(density_vector_s2); title=L"s_2"))    
end

if abspath(PROGRAM_FILE) != @__FILE__
    chain = h5open("Res/hpc-validation_chain-8.h5", "r") do f
        read(f, Chains)
    end

    a=draw_population_posterior(chain)
    # savefig(a, "d1_popu.eps")
end