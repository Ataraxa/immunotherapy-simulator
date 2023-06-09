% test
clc

params.td = 2.5;

params.k1 = 0.05;
params.k2 = 2.17;
params.k3 = 20;
params.k4 = 140;
params.k5 = 1.67;
params.k6 = 0.70;

params.d1 = 1.63;
params.d2 = 1.80;
params.d3 = 0.22;
params.d4 = 1.10;
params.d5 = 0.01;
params.d6 = 0.04;
params.d7 = 17;
params.d8 = 0.25;

params.s1 = 1.50;
params.s2 = 0.1;

save('struct_takuya.mat', 'params')
