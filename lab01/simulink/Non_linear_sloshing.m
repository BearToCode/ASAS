% Non linear pendolum
clc; clear;

addpath(genpath('.'));

params.density = 1000; % [kg/m^3]
params.h = 1; % [m]
params.d = 1; % [m]
params.g = 9.81; % [m/s^2]

damping = 0.003;
x0 = [ 0 ; 0 ]; % initial conditions

pendulums = sloshing_pendulums(params, 4);
%%
sim_nl = sim('Non_linear_sloshing.slx');

t = sim_nl.tout;     % Simulated time
u = sim_nl.u;        % Forcing term
theta= sim_nl.theta; % Problema: 3 colonne di theta sono vuote
force=sim_nl.force;
theta_deg=(theta./pi).*180;

figure
plot(t,force,'LineWidth',1);
xlabel('t [s]');
ylabel('$F_x$ [N]', 'Interpreter','LaTex')

figure
plot(t,theta_deg(:,1),t,theta_deg(:,2),t,theta_deg(:,3),t,theta_deg(:,4));
legend('I mode',"II mode","III mode","IV mode")
xlabel('t [s]');
ylabel('$\theta$ [°]', 'Interpreter','LaTex')
