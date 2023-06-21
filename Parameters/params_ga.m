    params_list = ["td", "t_delay", "t_last", "t_delay12", "t_last12", ...
        "k1", "k2", "k3", "k4", "k5", "k6", ...
        "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", ...
        "s1", "s2", ...
        "ig0", "c0", "p0", "vl0"];

i = 1;
params = struct;
for parameter = params1
    params.(params_list(i)) = parameter;
    i = i + 1;
end
disp(params)

save ("Parameters/struct_ga.mat", "params")
fprintf("Successfully saved new settings!\n")