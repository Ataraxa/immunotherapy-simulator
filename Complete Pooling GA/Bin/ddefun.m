function dydt= ddefun(t, y, Z)

% test
clc

t_in = 7;
t_delay = 0.587;
t_last = 3.5;

k1 = 0.05;
k2 = 2.17;
k3 = 83.9;
k4 = 1292;
k5 = 1.67;
k6 = 0.70;

d1 = 1.63;
d2 = 1.80;
d3 = 0.22;
d4 = 1.10;
d5 = 0.01;
d6 = 0.04;
d7 = 17;
d8 = 0.25;

s1 = 1.50;
s2 = 0.1;

v_max = 600;

if (t_in + t_delay < t) && (t < t_in + t_delay + t_last)
    d_cbd = 1;
else
    d_cbd = 0;
end

%MODEL_ODE Summary of this function goes here
%   Detailed explanation goes here
% TODO: add delay again
% g_lagged = Z(:,1); % what is going on here?
y_lagged = Z(:,1);

dydt = [
  k1 + k2*(d_cbd + 0) - d1*y(1); % dg/dt
  k3 + k4*y_lagged(1) - d2*y(2); % dc/dt
  k5 - (d3+d4*y(1))*y(3); % dp/dt
  k6*(1-(y(4)+y(5))/v_max)*y(4) - (d5 + (d6*y(2)/(1+s1*y(3)*(1-0)) + d7*y(1))/(1+s2*(y(4)+y(5))))*y(4); % dv_l/dt
  (d5 + (d6*y(2)/(1+s1*y(3)*(1-0)) + d7*y(1))/(1+s2*(y(4)+y(5))))*y(4) - d8*y(5);
];
end


