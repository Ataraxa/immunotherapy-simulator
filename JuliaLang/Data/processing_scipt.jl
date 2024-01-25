using DelimitedFiles

# Settings and data
data = readdlm("./Data/raw/full_data_day7.csv", ',')
n_mice = 5

function filter_matrix_col(col_index::Vector{Int64})
    return data[:,col_index]
end

function process_matrix(input_matrix)
    # Initialise final matrix
    output_matrix = zeros(Float64, n_mice, 7, 3)

    row_number = 1

    # Counters
    current_day = data[row_number, 2]
    forward_row = 0
    time = 1
    
    while time < 8
        while input_matrix[row_number + forward_row] == current_day
            output_matrix[forward_row+1,time,:] = input_matrix[row_number+forward_row,:]
        end
        time += 1
    end
end

middle_mat = filter_matrix_col([2, 5, 13, 16])