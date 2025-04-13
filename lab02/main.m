clc; clear;

addpath(genpath('.'));

params.M = 3; % [kg]
params.m = 0.1; % [kg]
params.l = 0.75; % [m]
params.g = 9.8; % [m/s^2]
params.c = 0.1; % [N/m/s]
params.b = 0.001; % [Nm/rad/s]

% Utility function for using arrayfun with cell arrays
map = @(f, x) cell2mat(arrayfun(@(idx) f(x(idx)), 1:length(x), "UniformOutput", false));

%% Task 1.1 – Mathematical model
% A. Derive the set of second-order nonlinear EOM of the
% system (show and comment every step of the derivation).
% B: B. Derive the nonlinear state-space formulation (show
% and comment every step of the derivation) with the force
% f as input and the cart position and pendulum angle as
% outputs.

[f, g] = ipend_nonlinear(params);

%% Task 1.2 – Nonlinear response to a force applied on the cart
% Starting from the initial position corresponding to the point x0 = [0; 0; pi; 0]T,
% develop a MATLAB code to simulate the nonlinear (open-loop) response in terms of cart
% position and pendulum angle in the time interval [0 – 20]s for the following three input cases:
% a) 1-second pulse force (from 1 to 2 s) on the cart of amplitude 1 N
% b) 1-second pulse force (from 1 to 2 s) on the cart of amplitude 5N
% c) 1-second pulse force (from 1 to 2 s) on the cart of amplitude 25 N
% A. Justify the choice of the integration method, by also showing and commenting its
% performance.
% B. Plot the output responses ( x(t) and theta(t) ) for each input case.
% C. Provide comments on the evolution of the variables for each input case and on the
% comparison of the responses in the three cases.

x0 = [0; 0; pi; 0]; % initial state [x; dx; theta; dtheta]
tspan = [0 20]; % time span for simulation

u1 = @(t) 1 * (t >= 1 & t <= 2); % 1 N pulse force
u2 = @(t) 5 * (t >= 1 & t <= 2); % 5 N pulse force
u3 = @(t) 25 * (t >= 1 & t <= 2); % 25 N pulse force

odefun1 = @(t, x) f(x, u1(t));
odefun2 = @(t, x) f(x, u2(t));
odefun3 = @(t, x) f(x, u3(t));

[t1, x1] = ode45(odefun1, tspan, x0);
[t2, x2] = ode45(odefun2, tspan, x0);
[t3, x3] = ode45(odefun3, tspan, x0);

y1 = map(@(idx) g(x1(idx, :), u1(t1(idx))), 1:length(t1));
y2 = map(@(idx) g(x2(idx, :), u2(t2(idx))), 1:length(t2));
y3 = map(@(idx) g(x3(idx, :), u3(t3(idx))), 1:length(t3));

pos1 = y1(1, :)'; % cart positions
pos2 = y2(1, :)';
pos3 = y3(1, :)';
theta1 = y1(2, :)' .* 180 / pi; % pendulum angles
theta2 = y2(2, :)' .* 180 / pi;
theta3 = y3(2, :)' .* 180 / pi;

task1_u1_pos = figure;
plot(t1, pos1, 'Color', "#77AC30");
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
legend('$x(t)$', 'Interpreter', 'latex');

save_figure(task1_u1_pos, 'task1_u1_pos.png')
title('Response to 1 N pulse force on the cart');

task1_u1_theta = figure;
plot(t1, theta1, 'Color', "#4DBEEE");
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
legend('$\theta(t)$', 'Interpreter', 'latex');

save_figure(task1_u1_theta, 'task1_u1_theta.png')
title('Response to 1 N pulse force on the cart');

task_u2_pos = figure;
plot(t2, pos2, 'Color', "#77AC30");
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
legend('$x(t)$', 'Interpreter', 'latex');

save_figure(task_u2_pos, 'task1_u2_pos.png')
title('Response to 5 N pulse force on the cart');

task_u2_theta = figure;
plot(t2, theta2, 'Color', "#4DBEEE");
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
legend('$\theta(t)$', 'Interpreter', 'latex');

save_figure(task_u2_theta, 'task1_u2_theta.png')
title('Response to 5 N pulse force on the cart');

task_u3_pos = figure;
plot(t3, pos3, 'Color', "#77AC30");
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
legend('$x(t)$', 'Interpreter', 'latex');

save_figure(task_u3_pos, 'task1_u3_pos.png')
title('Response to 25 N pulse force on the cart');

task_u3_theta = figure;
plot(t3, theta3, 'Color', "#4DBEEE");
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
legend('$\theta(t)$', 'Interpreter', 'latex');

save_figure(task_u3_theta, 'task1_u3_theta.png')
title('Response to 25 N pulse force on the cart');

%% Task 1.3 – Nonlinear response to a force applied on the cart
% A. Repeat Task 1.2 using a Simulink model based on integral block approach
