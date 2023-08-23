% Code to optimise the parameter vector
% Just needs the fitness function
clc;

% Load lower and upper bound
width_factor = 0.15;
load Parameters/struct_christian.mat params;
lower_bound = (1 - width_factor) .* struct2array(params);
upper_bound = (1 + width_factor) .* struct2array(params);

% Pre-GA debug
% empty

% Setting-up and using the GA built-in function
options = optimoptions('ga','UseParallel',true);
options.CrossoverFraction = 0.8;
opt_params = ga(@fitness_modCost, 25,...
    [], [], [], [],...
    lower_bound, ...  % lower bound for parameters optimisation
    upper_bound, ...  % upper bound for same
    [], options);

% opt_params = ga(@fitness_function_neutral, 25, [],[],[],[],[],[],[],options);
params_list = ["td", "t_delay", "t_last", "t_delay12", "t_last12", ...
        "k1", "k2", "k3", "k4", "k5", "k6", ...
        "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", ...
        "s1", "s2", ...
        "ig0", "c0", "p0", "vl0"];

i = 1;
clear("params")
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