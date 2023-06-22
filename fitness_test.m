clc; close all;

% load Parameters\struct_takuya.mat params;
load Archives\struct_ga_attempt1.mat params
params.ig0 = 0.015;
params.c0 = 10;
params.p0 = 4.4;
params.vl0 = 7.06;
disp(params)
vectorised_params = cell2mat(struct2cell(params));
tic;
fitness_function_neutral(vectorised_params);
toc
% set(gcf, 'Position', get(0, 'Screensize'));
fprintf("Test run successfully\n")