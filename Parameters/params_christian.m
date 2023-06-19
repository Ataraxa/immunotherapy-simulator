% test
clc; clear

params.td = 1.47;

params.t_delay = 0.38;
params.t_last = 5.41;

params.t_delay12 = 3.74;
params.t_last12 = 1.36;

params.k1 = 0.17;
params.k2 = 5.3;
params.k3 = 83.9;
params.k4 = 1292.16;
params.k5 = 4.85;
params.k6 = 0.49;

params.d1 = 9.97;
params.d2 = 11.05;
params.d3 = 1.09;
params.d4 = 5.8;
params.d5 = 0.02;
params.d6 = 0.04;
params.d7 = 51.94;
params.d8 = 0.61;

params.s1 = 16;
params.s2 = 0.35;
save('Parameters/struct_christian.mat', 'params')
fprintf("Successfully saved new settings!\n")
