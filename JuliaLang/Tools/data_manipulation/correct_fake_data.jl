using DelimitedFiles: readdlm, writedlm

data = readdlm("Data/fakeData/raw/trajectories-1.csv", ',')


function snapAtThreshold(x, thresh)
    return (x < thresh) ? thresh : x
end

data = snapAtThreshold.(data, 1e-6)

writedlm("Data/fakeData/trajectories-1.csv", data, ',')


