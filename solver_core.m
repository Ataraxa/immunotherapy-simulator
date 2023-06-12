% Solves the current model (modular parameters loading)
clc; clear; close all;

%% Settings
load("Parameters\struct_christian.mat"); % Select parameter set here
treatment = treatment_factory.cbd_9_14(); % Select treatment here

%% Solve and plot
sol = immuno_solver(params, treatment);

% Plotting results of the simulation
figure(1); hold on;
plot(sol.x,sol.y(4,:),'-','LineWidth',1)
plot(sol.x,sol.y(5,:),'-', 'LineWidth',1)
plot(sol.x,(sol.y(4,:)+sol.y(5,:)),'-','LineWidth',1)
xlabel('Day since tumour inoculation', 'Interpreter','latex');
ylabel('Tumour Volume ($mm^3$)', 'Interpreter','latex');
legend('Living Tumour', 'Dead Tumour','Total Tumour','Location','NorthWest');