params.td = 2.5;

params.t_delay = 0.57;
params.t_last = 3.80;
params.t_delay12 = 0.90;
params.t_last12 = 1.10;

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

params.ig0 = 0.015;
params.c0 = 12;
params.p0 = 4.4;
params.vl0 = 7.12;

save('Parameters/struct_takuya.mat', 'params')
