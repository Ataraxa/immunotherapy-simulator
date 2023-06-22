% Code to optimise the parameter vector
% Just needs the fitness function
clc;
options = optimoptions('ga','UseParallel',true);
opt_params = ga(@fitness_function_neutral, 25,...
    [], [], [], [],...
    [0.1,0.1,0.1,0.1,0.1, 0.01,0.1,20,500,0.1,0.1, 1,5,0.01,0.1,0.1,0.1,0.1,0.1, 0.1,0.1, 0.001,1,1,5], ...  % lower bound
    [10,10,10,10,10, 10,10,500,1500,10,10, 20,20,10,10,10,10,100,10, 100,10, 1,20,20,10], ...
    [], options);
params_list = ["td", "t_delay", "t_last", "t_delay12", "t_last12", ...
        "k1", "k2", "k3", "k4", "k5", "k6", ...
        "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", ...
        "s1", "s2", ...
        "ig0", "c0", "p0", "vl0"];

i = 1;
params = struct;    
for parameter = opt_params
    params.(params_list(i)) = parameter;
    i = i + 1;
end
disp(params)

save ("Parameters/struct_ga.mat", "params")
fprintf("Successfully saved new settings!\n")
% disp("starting top optimse")
% paramws = ga(@quick_testing, 21)

% [10,10,10,10,10, 10,10,500,1500,10,10, 20,20,10,10,10,10,100,10, 100,10, 1,20,20,10]