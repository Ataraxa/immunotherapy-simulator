clc; clear; 

load Parameters\struct_christian.mat params;
vectorised_params = cell2mat(struct2cell(params));
disp(vectorised_params)
fitness_function(vectorised_params)
fprintf("Test run successfully\n")