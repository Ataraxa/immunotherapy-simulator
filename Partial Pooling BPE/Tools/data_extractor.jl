using DelimitedFiles: readdlm

function read_data(selected_days, num_experiments, step_size)
    data_matrix = Array{Float64}(undef, num_experiments, length(selected_days))

    for i in 1:num_experiments
        data = readdlm("Data/trajectories-$i.csv", ',')
        exact_vol = data[5, :] + data[6, :]
        approx_vol = Array(exact_vol) + 10.0 * randn(size(exact_vol)[1])
        data_matrix[i,:] = approx_vol[selected_days .* trunc(Int, 1/step_size) .+ 1]
    end

    return(data_matrix)
end