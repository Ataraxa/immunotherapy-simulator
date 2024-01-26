using Distributions

christian_true_params = [
    1.872824160505723, 
    0.487923336960266,
    4.893475496158970, 
    3.697285769720824,
    1.075411285414352,
    # ⬇ index 6
    0.221936172352179,  #k1
    6.081963924603686,  #k2
    74.60149998881774,  #k3
    929.5659753695654,  #k4
    5.815052050675676,  #k5
    0.539494697127057,  #k6
    # ⬇ index 12
    10.177238034144443, # d1 
    10.610875554141064, # d2 
    1.3262458309759877, # d3
    4.4860998573217170, # d4
    0.0160245808033653, # d5
    0.0341621446708329, # d6
    59.611074450040185, # d7
    0.5697851077620799, # d8
    # ⬇ index 20
    14.068450976326906, # s1
    0.4108126097661595, # s2
    # ⬇ index 22
    0.007911968983770, # g (IFNγ)
    8.851474387488993, # c (CD8+)
    5.999229998549251, # p (PD-1)
    5.573047884761726  # v (Living tumour)
    # ⟹ Total: 25 params
]

function gen_priors(
    distro,
    std::Float64,
    is_info::Bool;
    base::Vector{Float64}=christian_true_params
    )

    return [truncated(distro((is_info ? log(par) : 0), std); lower=-7, upper=7)
         for par in base]
end

priors = gen_priors(Cauchy,1.,false)