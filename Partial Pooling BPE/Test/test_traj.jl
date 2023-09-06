#=
File to visualise the trajectory files and test their formatting options.
=#
import CSV, DelimitedFiles
using Plots: plot, plot!

function visualise_traj(filename)
    data = readdlm(filename,',')
    t = data[1,:]

    # plot(t, data[2,:])
    # plot!(t, data[3,:])
    # plot!(t, data[4,:])
    plot(t, data[5,:]+data[6,:])
end

visualise_traj("Data/trajectories-average.csv")