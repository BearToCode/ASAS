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

max_pos = max([max(pos1), max(pos2), max(pos3)]);
min_theta = min([min(theta1), min(theta2), min(theta3)]);
max_theta = max([max(theta1), max(theta2), max(theta3)]);

task1_u1_pos = figure;
plot(t1, pos1, 'Color', "#77AC30");
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos]);
legend('$x(t)$', 'Interpreter', 'latex');

save_figure(task1_u1_pos, 'task1_u1_pos.png')
title('Response to 1 N pulse force on the cart');

task1_u1_theta = figure;
plot(t1, theta1, 'Color', "#4DBEEE");
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta max_theta]);
legend('$\theta(t)$', 'Interpreter', 'latex');

save_figure(task1_u1_theta, 'task1_u1_theta.png')
title('Response to 1 N pulse force on the cart');

task_u2_pos = figure;
plot(t2, pos2, 'Color', "#77AC30");
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos]);
legend('$x(t)$', 'Interpreter', 'latex');

save_figure(task_u2_pos, 'task1_u2_pos.png')
title('Response to 5 N pulse force on the cart');

task_u2_theta = figure;
plot(t2, theta2, 'Color', "#4DBEEE");
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta max_theta]);
legend('$\theta(t)$', 'Interpreter', 'latex');

save_figure(task_u2_theta, 'task1_u2_theta.png')
title('Response to 5 N pulse force on the cart');

task_u3_pos = figure;
plot(t3, pos3, 'Color', "#77AC30");
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos]);
legend('$x(t)$', 'Interpreter', 'latex');

save_figure(task_u3_pos, 'task1_u3_pos.png')
title('Response to 25 N pulse force on the cart');

task_u3_theta = figure;
plot(t3, theta3, 'Color', "#4DBEEE");
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta max_theta]);
legend('$\theta(t)$', 'Interpreter', 'latex');

save_figure(task_u3_theta, 'task1_u3_theta.png')
title('Response to 25 N pulse force on the cart');

%% Task 1.3 – Nonlinear response to a force applied on the cart
% A. Repeat Task 1.2 using a Simulink model based on integral block approach
% B. Repeat Task 1.2 using a Simulink model based on M-function.

% TODO: add plots

%% Task 1.4 –Analysis of the linearized system
% A. Derive the mathematical expression and numerical values of the
% state-space matrices (A,B,C,D) of the linearized model about the
% equilibrium point when the outputs are the cart position and the pendulum angle.
% B. Derive the mathematical and numerical expression of the transfer functions (TFs).
% C. Plot the pole-zero map for both TFs.
% D. Plot the asymptotic and real Bode diagram (magnitude and phase) for both TFs.

linear_sys = ipend_linear(params);
[G_x, G_theta] = ipend_tf(params);

G_x_figure = figure;
pzmap(G_x);
save_figure(G_x_figure, 'task4_pzmap_x.png')
title('Pole-Zero Map of the Linearized System - $x$', 'Interpreter', 'latex');

G_theta_figure = figure;
pzmap(G_theta);
save_figure(G_theta_figure, 'task4_pzmap_theta.png')
title('Pole-Zero Map of the Linearized System - $\theta$', 'Interpreter', 'latex');

db = @(x) 20 * log10(x);

w = logspace(-3, 2, 1000); % frequency range for Bode plot
[mag_asymp_x, phase_asymp_x] = asymp_bode(G_x, w);
[mag_asymp_theta, phase_asymp_theta] = asymp_bode(G_theta, w);

[mag_real_x, phase_real_x] = bode(G_x, w);
[mag_real_theta, phase_real_theta] = bode(G_theta, w);

% Convert the three dimensional arrays to one dimension
mag_real_x = squeeze(mag_real_x);
mag_real_theta = squeeze(mag_real_theta);
phase_real_x = squeeze(phase_real_x);
phase_real_theta = squeeze(phase_real_theta);

task4_bode_x = figure;
subplot(2, 1, 1);
semilogx(w, db(mag_asymp_x), 'DisplayName', 'Asymptotic');
hold on;
semilogx(w, db(mag_real_x), 'DisplayName', 'Real');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
legend('$|\bar{G_x}(j\omega)|$', '$|G_x(j\omega)|$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_x$', 'Interpreter', 'latex');

subplot(2, 1, 2);
semilogx(w, phase_asymp_x, 'DisplayName', 'Asymptotic');
hold on;
semilogx(w, phase_real_x, 'DisplayName', 'Real');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [°]');
legend('$\angle \bar{G_x}(j\omega)$', '$\angle G_x(j\omega)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_x$', 'Interpreter', 'latex');

save_figure(task4_bode_x, 'task4_bode_x.png')

task4_bode_theta = figure;
subplot(2, 1, 1);
semilogx(w, db(mag_asymp_theta), 'DisplayName', 'Asymptotic');
hold on;
semilogx(w, db(mag_real_theta), 'DisplayName', 'Real');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
legend('$|\bar{G_\theta}(j\omega)|$', '$|G_\theta(j\omega)|$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_\theta$', 'Interpreter', 'latex');

subplot(2, 1, 2);
semilogx(w, phase_asymp_theta, 'DisplayName', 'Asymptotic');
hold on;
semilogx(w, phase_real_theta, 'DisplayName', 'Real');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [°]');
legend('$\angle \bar{G_\theta}(j\omega)$', '$\angle G_\theta(j\omega)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_\theta$', 'Interpreter', 'latex');

save_figure(task4_bode_theta, 'task4_bode_theta.png')

%% Task 1.5 – Comparison of the nonlinear and linear response
% Starting from the initial position corresponding to the equilibrium point x0
% develop a MATLAB code to simulate the linear (open-loop) response in terms of cart
% position and pendulum angle in the time interval [0 – 20] s for the following three input cases:
% a) 1-second pulse force (from 1 to 2 s) on the cart of amplitude 1 N
% b) 1-second pulse force (from 1 to 2 s) on the cart of amplitude 5N
% c) 1-second pulse force (from 1 to 2 s) on the cart of amplitude 25 N
% For each input case, plot on the same figure the linear response and the nonlinear response
% computed in Task 1.2 (x(t) and theta(t)) and provide comments on the differences.

odefun1 = @(t, x) linear_sys.A * x + linear_sys.B * u1(t);
odefun2 = @(t, x) linear_sys.A * x + linear_sys.B * u2(t);
odefun3 = @(t, x) linear_sys.A * x + linear_sys.B * u3(t);

[t1_lin, x1_lin] = ode45(odefun1, tspan, zeros(4, 1));
[t2_lin, x2_lin] = ode45(odefun2, tspan, zeros(4, 1));
[t3_lin, x3_lin] = ode45(odefun3, tspan, zeros(4, 1));

y1_lin = linear_sys.C * x1_lin' + linear_sys.D * u1(t1_lin)';
y2_lin = linear_sys.C * x2_lin' + linear_sys.D * u2(t2_lin)';
y3_lin = linear_sys.C * x3_lin' + linear_sys.D * u3(t3_lin)';

pos1_lin = y1_lin(1, :)'; % cart positions
pos2_lin = y2_lin(1, :)';
pos3_lin = y3_lin(1, :)';
theta1_lin = y1_lin(2, :)' .* 180 / pi + 180; % pendulum angles
theta2_lin = y2_lin(2, :)' .* 180 / pi + 180;
theta3_lin = y3_lin(2, :)' .* 180 / pi + 180;

max_pos_lin = max([max(pos1_lin), max(pos2_lin), max(pos3_lin)]);
min_theta_lin = min([min(theta1_lin), min(theta2_lin), min(theta3_lin)]);
max_theta_lin = max([max(theta1_lin), max(theta2_lin), max(theta3_lin)]);

task5_u1_pos_lin = figure;
plot(t1, pos1, 'Color', "#77AC30", 'DisplayName', 'Nonlinear');
hold on;
plot(t1_lin, pos1_lin, 'Color', "#D95319", 'DisplayName', 'Linear');
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos_lin]);
legend('$x(t)$', '$\bar{x}(t)$', 'Interpreter', 'latex');

save_figure(task5_u1_pos_lin, 'task5_u1_pos_lin.png')
title('Response to 1 N pulse force on the cart');

task5_u1_theta_lin = figure;
plot(t1, theta1, 'Color', "#4DBEEE", 'DisplayName', 'Nonlinear');
hold on;
plot(t1_lin, theta1_lin, 'Color', "#A2142F", 'DisplayName', 'Linear');
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta_lin max_theta_lin]);
legend('$\theta(t)$', '$\bar{\theta}(t)$', 'Interpreter', 'latex');

save_figure(task5_u1_theta_lin, 'task5_u1_theta_lin.png')
title('Response to 1 N pulse force on the cart');

task5_u2_pos_lin = figure;
plot(t2, pos2, 'Color', "#77AC30", 'DisplayName', 'Nonlinear');
hold on;
plot(t2_lin, pos2_lin, 'Color', "#D95319", 'DisplayName', 'Linear');
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos_lin]);
legend('$x(t)$', '$\bar{x}(t)$', 'Interpreter', 'latex');

save_figure(task5_u2_pos_lin, 'task5_u2_pos_lin.png')
title('Response to 5 N pulse force on the cart');

task5_u2_theta_lin = figure;
plot(t2, theta2, 'Color', "#4DBEEE", 'DisplayName', 'Nonlinear');
hold on;
plot(t2_lin, theta2_lin, 'Color', "#A2142F", 'DisplayName', 'Linear');
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta_lin max_theta_lin]);
legend('$\theta(t)$', '$\bar{\theta}(t)$', 'Interpreter', 'latex');

save_figure(task5_u2_theta_lin, 'task5_u2_theta_lin.png')
title('Response to 5 N pulse force on the cart');

task5_u3_pos_lin = figure;
plot(t3, pos3, 'Color', "#77AC30", 'DisplayName', 'Nonlinear');
hold on;
plot(t3_lin, pos3_lin, 'Color', "#D95319", 'DisplayName', 'Linear');
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos_lin]);
legend('$x(t)$', '$\bar{x}(t)$', 'Interpreter', 'latex');

save_figure(task5_u3_pos_lin, 'task5_u3_pos_lin.png')
title('Response to 25 N pulse force on the cart');

task5_u3_theta_lin = figure;
plot(t3, theta3, 'Color', "#4DBEEE", 'DisplayName', 'Nonlinear');
hold on;
plot(t3_lin, theta3_lin, 'Color', "#A2142F", 'DisplayName', 'Linear');
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta_lin max_theta_lin]);
legend('$\theta(t)$', '$\bar{\theta}(t)$', 'Interpreter', 'latex');

save_figure(task5_u3_theta_lin, 'task5_u3_theta_lin.png')
title('Response to 25 N pulse force on the cart');

%% Task 2.1 – Mathematical model
% A. Derive the nonlinear EOM of the system
% (show and comment every step of the derivation).
% B. Derive the nonlinear state-space formulation with the force u
% as control input, the force d as disturbance input, and the cart
% position and pendulum angle as outputs (show and comment every
% step of the derivation).
% C. Derive the linearized EOM of the system (show and comment
% every step of the derivation).
% D. Derive the linearized state-space formulation (show and comment
% every step of the derivation).

clear; % clear workspace

params.M = 3; % [kg]
params.m = 0.1; % [kg]
params.l = 0.75; % [m]
params.g = 9.8; % [m/s^2]

[A, B_u, B_d, C, D_u, D_d] = ipend_control_linear(params);

%% Task 2.2 – Open-loop transfer function
% A. Derive the mathematical and numerical expression of the open-loop
% transfer function G(s) from the control input to the pendulum angle
% B. Derive the mathematical expression and compute the numerical values of
% the zeros and poles of G(s) and plot the pole/zero map.

[~, G] = ipend_control_tf(params);

G_figure = figure;
pzmap(G);
save_figure(G_figure, 'task2_pzmap.png')
title('Pole-Zero Map of the Open-Loop Transfer Function - $G(s)$', 'Interpreter', 'latex');

%% Task 2.3 – Proportional (P) control for stabilization
% Let’s consider a proportional (P) control:
% A. By assuming a reference equal to zero (r(s) = 0), draw the block diagram of the closed-loop system and
% derive the mathematical expression of the closed-loop transfer function from the disturbance input to
% the pendulum angle output
% B. Demonstrate analytically that a P control cannot stabilize the system.
% C. Confirm the previous result using the root locus technique and discuss the loci behaviour.
% D. Determine the minimum value of the control gain such that the closed-loop system is marginally stable.
