clc; clear; close all;

load Parameters\struct_christian.mat params;
vectorised_params = cell2mat(struct2cell(params));
fitness_function(vectorised_params)
set(gcf, 'Position', get(0, 'Screensize'));
fprintf("Test run successfully\n")