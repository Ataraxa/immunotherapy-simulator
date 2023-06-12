% Code to optimise the parameter vector
% Just needs the fitness function
clc;
params = ga(@fitness_function, 21,...
    [], [], [], [],...
    [0.1,0.1,0.1,0.1,0.1, 0,0,20,500,0,0, 0,0,0,0,0,0,0,0, 0,0,0 ],...  % lower bound
    [10,10,10,10,10, 10,10,500,1500,10,10, 20,20,10,10,10,10,100,10, 20,10,20]); % upper bound
% disp("starting top optimse")
% paramws = ga(@quick_testing, 21)