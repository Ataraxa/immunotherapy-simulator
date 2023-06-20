clear; close all;

load Parameters\struct_christian.mat params;
vectorised_params = cell2mat(struct2cell(params));
tic;
fitness_function(vectorised_params);
toc
% set(gcf, 'Position', get(0, 'Screensize'));
fprintf("Test run successfully\n")