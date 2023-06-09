% Solves the current model (modular parameters loading)
clc; clear; close all;
load("Parameters\struct_christian.mat");
treatment = treatment_factory.il_12();

sol = param_solver(params, treatment);

% Plotting results of the simulation
figure(1); hold on;
plot(sol.x,sol.y(4,:),'-')
plot(sol.x,sol.y(5,:),'-')
plot(sol.x,(sol.y(4,:)+sol.y(5,:)),'-')
xlabel('Day since tumour inoculation', 'Interpreter','latex');
ylabel('Tumour Volume ($mm^3$)', 'Interpreter','latex');
legend('Living Tumour', 'Dead Tumour','Total Tumour','Location','NorthWest');