using DelimitedFiles
using Plots: plot, scatter

data = readdlm("Data/experimental/tumour_vol_data.csv", ',')
processed = Dict()

### Settings
max_days = 30; max_nb_patiens = 30

### Process each row individually 
for data_point in 2:size(data)[1]
    tr_name = data[data_point, 1]
    day = data[data_point, 2]
    vol = data[data_point, 3]

    # Initialise matrix if first encounter 
    if !haskey(processed, tr_name)
        processed[tr_name] = Dict()
        processed[tr_name]["days"] = [day]
        processed[tr_name]["vol"] = zeros(max_days, max_nb_patiens) # row x col 
    end

    # Update volume info 
    for patient in 1:max_nb_patiens
        if  processed[tr_name]["vol"][patient, day] == 0
            processed[tr_name]["vol"][patient, day] = vol 
            break
        end
    end

    # Update valid days array 
    if !(day in processed[tr_name]["days"])
        push!(processed[tr_name]["days"], day)
    end
end

### Reprocess entries and plot
for treatment in keys(processed)
    slice = processed[treatment]["vol"]
    slice = slice[:, vec(mapslices(col -> any(col .!= 0), slice, dims = 1))]
    slice = slice[vec(mapslices(col -> any(col .!= 0), slice, dims = 2)), :]    
    display(plot(processed[treatment]["days"], slice'; title=treatment, legend=true))
end

