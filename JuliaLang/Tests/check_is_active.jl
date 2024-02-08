using Plots
include("../Model/mechanistic_model.jl")
include("../Model/treatments_lib.jl")

t_array = 0:0.1:27.0
res_cbd = zeros(length(t_array))
res_il12 = zeros(length(t_array))
tr = CBD_IL_12_ver7
t_delay = 0.487923336960266
t_delay12 = 1.075411285414352
t_last = 4.893475496158970
t_last12 = 1.075411285414352
for (i,t) in enumerate(t_array)
    res_cbd[i] = check_active(t, tr.t_in, t_delay, t_last, (tr.t_in != 0))
    res_il12[i] = check_active(t, tr.t_in12, t_delay12, t_last12, (tr.t_in12 != 0))
end
plt1 = plot()
plot!(t_array, res_il12)
display(plt1)
