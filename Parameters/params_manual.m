params.td = 8.3;
params.t_delay = 0.3;
params.t_last = 3.7;
params.t_delay12 = 9;
params.t_last12 = 2.1;
params.k1 = 0.01;
params.k2 = 3.1;
params.k3 = 20.5;
params.k4 = 1474.8;
params.k5 = 0.4;
params.k6 = 6.8;
params.d1 = 4.3;
params.d2 = 6;
params.d3 = 1.3;
params.d4 = 4.5;
params.d5 = 4.8;
params.d6 = 4.3;
params.d7 = 2.8;
params.d8 = 2.4;
params.s1 = 0.4;
params.s2 = 6.2;

params.ig0 = 0.001;
params.c0 = 3;
params.p0 = 4.4;
params.vl0 = 5;

save('Parameters/struct_manual.mat', 'params')
fprintf("Successfully saved new settings!\n")
