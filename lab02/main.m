clc; clear;

addpath(genpath('.'));
addpath(genpath('../lib'));

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

u1_color = "#0c66e4";
u1_color_accent = "#579dff";
u2_color = "#1f845a";
u2_color_accent = "#4bce97";
u3_color = "#c9372c";
u3_color_accent = "#f87168";

figure;
subplot(2, 3, 1)
plot(t1, pos1, 'Color', u1_color, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos]);
legend('$x(t)$', 'Interpreter', 'latex');
title('$x(t)$ response to $u_1(t)$', 'Interpreter', 'latex');

subplot(2, 3, 2)
plot(t2, pos2, 'Color', u2_color, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos]);
legend('$x(t)$', 'Interpreter', 'latex');
title('$x(t)$ response to $u_2(t)$', 'Interpreter', 'latex');

subplot(2, 3, 3)
plot(t3, pos3, 'Color', u3_color, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos]);
legend('$x(t)$', 'Interpreter', 'latex');
title('$x(t)$ response to $u_3(t)$', 'Interpreter', 'latex');

subplot(2, 3, 4)
plot(t1, theta1, 'Color', u1_color, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta max_theta]);
legend('$\theta(t)$', 'Interpreter', 'latex');
title('$\theta(t)$ response to $u_1(t)$', 'Interpreter', 'latex');

subplot(2, 3, 5)
plot(t2, theta2, 'Color', u2_color, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta max_theta]);
legend('$\theta(t)$', 'Interpreter', 'latex');
title('$\theta(t)$ response to $u_2(t)$', 'Interpreter', 'latex');

subplot(2, 3, 6)
plot(t3, theta3, 'Color', u3_color, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta max_theta]);
legend('$\theta(t)$', 'Interpreter', 'latex');
title('$\theta(t)$ response to $u_3(t)$', 'Interpreter', 'latex');

save_figure('task1_composite.png', true)

% Study the numerical performance of the integration method
t1_ode45 = t1; t2_ode45 = t2; t3_ode45 = t3;
pos1_ode45 = pos1; pos2_ode45 = pos2; pos3_ode45 = pos3;
theta1_ode45 = theta1; theta2_ode45 = theta2; theta3_ode45 = theta3;

% Use ode89 as reference method to estimate error
[t1_ode89, x1_ode89] = ode89(odefun1, tspan, x0); % using ode89
[t2_ode89, x2_ode89] = ode89(odefun2, tspan, x0);
[t3_ode89, x3_ode89] = ode89(odefun3, tspan, x0);

pos1_ode89 = x1_ode89(:, 1); % cart positions
pos2_ode89 = x2_ode89(:, 1);
pos3_ode89 = x3_ode89(:, 1);

theta1_ode89 = x1_ode89(:, 3) .* 180 / pi; % pendulum angles
theta2_ode89 = x2_ode89(:, 3) .* 180 / pi;
theta3_ode89 = x3_ode89(:, 3) .* 180 / pi;

dt1_ode45 = diff(t1_ode45); % time step for ode45
dt2_ode45 = diff(t2_ode45);
dt3_ode45 = diff(t3_ode45);

dt1_ode45_mov50 = movmean(dt1_ode45, 50); % moving average of the time step for ode45
dt2_ode45_mov50 = movmean(dt2_ode45, 50);
dt3_ode45_mov50 = movmean(dt3_ode45, 50);

f_error = @(t, value, ref_t, ref_values) abs(value - interp1(ref_t, ref_values, t)); % error function

pos1_error = f_error(t1_ode45, pos1_ode45, t1_ode89, pos1_ode89); % error for cart position
pos2_error = f_error(t2_ode45, pos2_ode45, t2_ode89, pos2_ode89);
pos3_error = f_error(t3_ode45, pos3_ode45, t3_ode89, pos3_ode89);

theta1_error = f_error(t1_ode45, theta1_ode45, t1_ode89, theta1_ode89); % error for pendulum angle
theta2_error = f_error(t2_ode45, theta2_ode45, t2_ode89, theta2_ode89);
theta3_error = f_error(t3_ode45, theta3_ode45, t3_ode89, theta3_ode89);

max_dt = max([max(dt1_ode45), max(dt2_ode45), max(dt3_ode45)]); % maximum time step for ode45
max_pos_error = max([max(pos1_error), max(pos2_error), max(pos3_error)]); % maximum error for cart position
max_theta_error = max([max(theta1_error), max(theta2_error), max(theta3_error)]); % maximum error for pendulum angle

% Create a 3x3 grid of subplots
figure;
subplot(3, 3, 1);
plot(dt1_ode45, 'Color', u1_color, 'DisplayName', 'dt (ode45)', 'LineWidth', 1.5);
hold on;
plot(dt1_ode45_mov50, 'Color', u1_color_accent, 'DisplayName', 'dt (ode45 movmean)', 'LineWidth', 1.5);
grid on;
xlabel('Index');
ylabel('Time step [s]');
ylim([0 max_dt]);
legend('dt (ode45)', 'dt (ode45 movmean)', 'Location', 'Best');
title('Time step for $u_1(t)$', 'Interpreter', 'latex');

subplot(3, 3, 2);
plot(dt2_ode45, 'Color', u2_color, 'DisplayName', 'dt (ode45)', 'LineWidth', 1.5);
hold on;
plot(dt2_ode45_mov50, 'Color', u2_color_accent, 'DisplayName', 'dt (ode45 movmean)', 'LineWidth', 1.5);
grid on;
xlabel('Index');
ylabel('Time step [s]');
ylim([0 max_dt]);
legend('dt (ode45)', 'dt (ode45 movmean)', 'Location', 'Best');
title('Time step for $u_2(t)$', 'Interpreter', 'latex');

subplot(3, 3, 3);
plot(dt3_ode45, 'Color', u3_color, 'DisplayName', 'dt (ode45)', 'LineWidth', 1.5);
hold on;
plot(dt3_ode45_mov50, 'Color', u3_color_accent, 'DisplayName', 'dt (ode45 movmean)', 'LineWidth', 1.5);
grid on;
xlabel('Index');
ylabel('Time step [s]');
ylim([0 max_dt]);
legend('dt (ode45)', 'dt (ode45 movmean)', 'Location', 'Best');
title('Time step for $u_3(t)$', 'Interpreter', 'latex');

subplot(3, 3, 4);
plot(t1_ode45, pos1_error, 'Color', u1_color, 'DisplayName', 'Error (ode45)', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error [m]');
ylim([0 max_pos_error]);
legend('Error (ode45)', 'Location', 'Best');
title('$x(t)$ error for $u_1(t)$', 'Interpreter', 'latex');

subplot(3, 3, 5);
plot(t2_ode45, pos2_error, 'Color', u2_color, 'DisplayName', 'Error (ode45)', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error [m]');
ylim([0 max_pos_error]);
legend('Error (ode45)', 'Location', 'Best');
title('$x(t)$ error for $u_2(t)$', 'Interpreter', 'latex');

subplot(3, 3, 6);
plot(t3_ode45, pos3_error, 'Color', u3_color, 'DisplayName', 'Error (ode45)', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error [m]');
ylim([0 max_pos_error]);
legend('Error (ode45)', 'Location', 'Best');
title('$x(t)$ error for $u_3(t)$', 'Interpreter', 'latex');

subplot(3, 3, 7);
plot(t1_ode45, theta1_error, 'Color', u1_color, 'DisplayName', 'Error (ode45)', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error [°]');
ylim([0 max_theta_error]);
legend('Error (ode45)', 'Location', 'Best');
title('$\theta(t)$ error for $u_1(t)$', 'Interpreter', 'latex');

subplot(3, 3, 8);
plot(t2_ode45, theta2_error, 'Color', u2_color, 'DisplayName', 'Error (ode45)', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error [°]');
ylim([0 max_theta_error]);
legend('Error (ode45)', 'Location', 'Best');
title('$\theta(t)$ error for $u_2(t)$', 'Interpreter', 'latex');

subplot(3, 3, 9);
plot(t3_ode45, theta3_error, 'Color', u3_color, 'DisplayName', 'Error (ode45)', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error [°]');
ylim([0 max_theta_error]);
legend('Error (ode45)', 'Location', 'Best');
title('$\theta(t)$ error for $u_3(t)$', 'Interpreter', 'latex');

save_figure('task1_errors.png', true)

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

figure;
subplot(1, 2, 1);
pzmap(G_x);
title('Pole-Zero Map of the Linearized System - $x$', 'Interpreter', 'latex');

subplot(1, 2, 2);
pzmap(G_theta);
title('Pole-Zero Map of the Linearized System - $\theta$', 'Interpreter', 'latex');

save_figure('task4_pzmap_composite.png', true)

Co = ctrb(linear_sys);

if det(Co) == 0
    disp("The system is not controllable")
end

Ob_x = obsv(linear_sys.A, linear_sys.C(1, :));

if det(Ob_x) == 0
    disp("x is not comp. observable")
end

Ob_t = obsv(linear_sys.A, linear_sys.C(2, :));

if det(Ob_t) == 0
    disp("theta is not comp. observable")
end

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

figure;
subplot(2, 1, 1);
semilogx(w, db(mag_asymp_x), 'DisplayName', 'Asymptotic', 'LineWidth', 1.5);
hold on;
semilogx(w, db(mag_real_x), 'DisplayName', 'Real', 'LineWidth', 1.5);
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
legend('$|\bar{G_x}(j\omega)|$', '$|G_x(j\omega)|$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_x$', 'Interpreter', 'latex');

subplot(2, 1, 2);
semilogx(w, phase_asymp_x, 'DisplayName', 'Asymptotic', 'LineWidth', 1.5);
hold on;
semilogx(w, phase_real_x, 'DisplayName', 'Real', 'LineWidth', 1.5);
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [°]');
legend('$\angle \bar{G_x}(j\omega)$', '$\angle G_x(j\omega)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_x$', 'Interpreter', 'latex');

save_figure('task4_bode_x.png')

figure;
subplot(2, 1, 1);
semilogx(w, db(mag_asymp_theta), 'DisplayName', 'Asymptotic', 'LineWidth', 1.5);
hold on;
semilogx(w, db(mag_real_theta), 'DisplayName', 'Real', 'LineWidth', 1.5);
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
legend('$|\bar{G_\theta}(j\omega)|$', '$|G_\theta(j\omega)|$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_\theta$', 'Interpreter', 'latex');

subplot(2, 1, 2);
semilogx(w, phase_asymp_theta, 'DisplayName', 'Asymptotic', 'LineWidth', 1.5);
hold on;
semilogx(w, phase_real_theta, 'DisplayName', 'Real', 'LineWidth', 1.5);
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [°]');
legend('$\angle \bar{G_\theta}(j\omega)$', '$\angle G_\theta(j\omega)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_\theta$', 'Interpreter', 'latex');

save_figure('task4_bode_theta.png')

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

figure;
subplot(2, 3, 1)
plot(t1, pos1, 'Color', u1_color, 'DisplayName', 'Nonlinear', 'LineWidth', 1.5);
hold on;
plot(t1_lin, pos1_lin, 'Color', u1_color_accent, 'DisplayName', 'Linear', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos_lin]);
legend('$x(t)$', '$\bar{x}(t)$', 'Interpreter', 'latex');
title('$x(t)$ linear response to $u_1(t)$', 'Interpreter', 'latex');

subplot(2, 3, 2)
plot(t2, pos2, 'Color', u2_color, 'DisplayName', 'Nonlinear', 'LineWidth', 1.5);
hold on;
plot(t2_lin, pos2_lin, 'Color', u2_color_accent, 'DisplayName', 'Linear', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos_lin]);
legend('$x(t)$', '$\bar{x}(t)$', 'Interpreter', 'latex');
title('$x(t)$ linear response to $u_2(t)$', 'Interpreter', 'latex');

subplot(2, 3, 3)
plot(t3, pos3, 'Color', u3_color, 'DisplayName', 'Nonlinear', 'LineWidth', 1.5);
hold on;
plot(t3_lin, pos3_lin, 'Color', u3_color_accent, 'DisplayName', 'Linear', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Cart position [m]');
ylim([0 max_pos_lin]);
legend('$x(t)$', '$\bar{x}(t)$', 'Interpreter', 'latex');
title('$x(t)$ linear response to $u_3(t)$', 'Interpreter', 'latex');

subplot(2, 3, 4)
plot(t1, theta1, 'Color', u1_color, 'DisplayName', 'Nonlinear', 'LineWidth', 1.5);
hold on;
plot(t1_lin, theta1_lin, 'Color', u1_color_accent, 'DisplayName', 'Linear', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta_lin max_theta_lin]);
legend('$\theta(t)$', '$\bar{\theta}(t)$', 'Interpreter', 'latex');
title('$\theta(t)$ linear response to $u_1(t)$', 'Interpreter', 'latex');

subplot(2, 3, 5)
plot(t2, theta2, 'Color', u2_color, 'DisplayName', 'Nonlinear', 'LineWidth', 1.5);
hold on;
plot(t2_lin, theta2_lin, 'Color', u2_color_accent, 'DisplayName', 'Linear', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta_lin max_theta_lin]);
legend('$\theta(t)$', '$\bar{\theta}(t)$', 'Interpreter', 'latex');
title('$\theta(t)$ linear response to $u_2(t)$', 'Interpreter', 'latex');

subplot(2, 3, 6)
plot(t3, theta3, 'Color', u3_color, 'DisplayName', 'Nonlinear', 'LineWidth', 1.5);
hold on;
plot(t3_lin, theta3_lin, 'Color', u3_color_accent, 'DisplayName', 'Linear', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Pendulum angle [°]');
ylim([min_theta_lin max_theta_lin]);
legend('$\theta(t)$', '$\bar{\theta}(t)$', 'Interpreter', 'latex');
title('$\theta(t)$ linear response to $u_3(t)$', 'Interpreter', 'latex');

save_figure('task5_composite_lin.png', true)

%%  Task 1.6 – Comparison of the nonlinear and linear response
% Repeat Task 1.5 using Simulink models

for k = [1, 5, 25]

    sim_pulse = sim('simulink_double_step.slx');

    u = sim_pulse.u; % Forcing term
    y = sim_pulse.y; % Output history

    hold on;
    subplot(1, 3, 1)
    plot(y.Time, u)
    xlabel("Time [s]")
    ylabel("Input [N]")
    hold on;
    subplot(1, 3, 2)
    plot(y.Time, y.Data(:, 1))
    xlabel("Time [s]")
    ylabel("x [m]")
    hold on;
    subplot(1, 3, 3)
    plot(y.Time, y.Data(:, 2))
    xlabel("Time [s]")
    ylabel("Theta [-]")

end

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

[f, g] = ipend_control_nonlinear(params);
[A, B_u, B_d, C, D_u, D_d] = ipend_control_linear(params);

%% Task 2.2 – Open-loop transfer function
% A. Derive the mathematical and numerical expression of the open-loop
% transfer function G(s) from the control input to the pendulum angle
% B. Derive the mathematical expression and compute the numerical values of
% the zeros and poles of G(s) and plot the pole/zero map.

[G_x, G_theta] = ipend_control_tf(params);

figure;
pzmap(G_theta);
title('Pole-Zero Map of the Open-Loop Transfer Function - $G(s)$', 'Interpreter', 'latex');
save_figure('task2_pzmap.png')

%% Task 2.3 – Proportional (P) control for stabilization
% Let’s consider a proportional (P) control:
% A. By assuming a reference equal to zero (r(s) = 0), draw the block diagram of the closed-loop system and
% derive the mathematical expression of the closed-loop transfer function from the disturbance input to
% the pendulum angle output
% B. Demonstrate analytically that a P control cannot stabilize the system.
% C. Confirm the previous result using the root locus technique and discuss the loci behavior.
% D. Determine the minimum value of the control gain such that the closed-loop system is marginally stable.

rlocus(G_theta);

% TODO: D

%% Task 2.4 – Proportional-Derivative (PD) control for stabilization
% Let’s consider a proportional-derivative (PD) control:
% A. By assuming a reference equal to zero (r(s) = 0), draw the block
% diagram of the closed-loop system and derive the mathematical
% expression of the closed-loop transfer function from the
% disturbance to the pendulum angle
% B. Determine the value of the control gains kp and kd such that
% the closed-loop pendulum response to a step disturbance
% would have peak-time less than 1 s and overshoot less than
% 20% (make a connection of these specifications to the
% locations of desired closed–loop poles)
% C. (MATLAB® code) Simulate and report the nonlinear
% closed-loop time response (cart position and pendulum
% angle) of the PD-controlled system with the gains
% obtained from point B., when the cart is subjected to a
% pulse disturbance force of 1 N and duration 0.1 s. Report
% also the time history of the control force and comment
% the behavior with respect to the system response.
% D. Repeat point C. using a Simulink® model.

Kp_theta = 606;
Kd_theta = 47;

% Integration parameters
x0 = [0; 0; 0; 0];
r = @(t) x0; % reference input
d = @(t) 1 * (t <= 0.1);
tspan = [0 20];

% Partial state PD controller
pd_theta = @(x) Kp_theta * x(3) + Kd_theta * x(4); % PD control law

odefun = @(t, x) f(x, pd_theta(r(t) - x), d(t));
[t, x] = ode45(odefun, tspan, x0);

theta = x(:, 3) .* 180 / pi;
control = arrayfun(@(idx) pd_theta(r(t(idx)) - x(idx, :)'), length(t)); % control force

max_theta = max(abs(theta));
max_control = max(abs(control));

figure;
plot(t, control, 'DisplayName', 'Control Force');
hold on;
grid on;
xlabel('Time [s]');
ylabel('Control Force [N]');
ylim([-max_control max_control]);
yyaxis right;
plot(t, theta, 'DisplayName', 'Pendulum Angle');
ylabel('Pendulum Angle [°]');
ylim([-max_theta max_theta]);
legend('$F_c(t)$', '$\theta(t)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Control Force and Pendulum Angle - PD Control', 'Interpreter', 'latex');
% MANCA IL PLOT DI X!

save_figure('task2_control_theta.png')

% D: Simulink

sim_control = sim('nonlinear_control.slx');

d_sim = sim_control.d; % Disturbance
x_sim = sim_control.x; % Output history
u_sim = sim_control.u; % Control force


figure
plot(x_sim.time, x_sim.signals.values(:, 1));
xlabel("Time [s]");
ylabel('Cart position [m]');

max_theta_sim = max(abs(x_sim.signals.values(:,3)*180/pi));
max_control_sim = max(abs(u_sim.signals.values));

figure;
plot(u_sim.time, u_sim.signals.values,'DisplayName', 'Control Force');
hold on;
grid on;
xlabel('Time [s]');
ylabel('Control Force [N]');
ylim([-max_control_sim max_control_sim]);
yyaxis right;
plot(x_sim.time,  x_sim.signals.values(:, 3) * 180 / pi, 'DisplayName', 'Pendulum Angle');
ylabel('Pendulum Angle [°]');
ylim([-max_theta_sim max_theta_sim]);
legend('$F_c(t)$', '$\theta(t)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Control Force and Pendulum Angle - PD Control', 'Interpreter', 'latex');


%% Task 2.5 – Partial and full state feedback control
% A. Show that the PD controller designed in Task 2.4 corresponds to a partial state feedback control, write
% the corresponding closed-loop linear state-space model, evaluate the resulting closed-loop poles
% and show the limitations of the PD control solution in terms of feasibility on a real system.
% B. Extending the PD approach to a full state feedback control, compute the control gains and evaluate
% the effect on performance (cart position and pendulum angle) and control effort of different selections
% of closed-loop poles. Comment the differences with respect to point A.

% Create a surface for the closed-loop system showing the max real part of the roots
Kp_x_space = linspace(-100, 100, 50);
Kd_x_space = linspace(-100, 100, 50);

[Kp_x, Kd_x] = meshgrid(Kp_x_space, Kd_x_space);
[F_den] = ipend_full_control_den(params);

f_plot = @(Kp_x, Kd_x) max(real(roots(F_den(Kp_x, Kd_x, Kp_theta, Kd_theta))));

z = arrayfun(f_plot, Kp_x, Kd_x);

figure;
surf(Kp_x, Kd_x, z);
xlabel('Kp_x');
ylabel('Kd_x');
zlabel('Real Part of Roots');
title('Real Part of Roots of the Closed-Loop System', 'Interpreter', 'latex');
view(235, 45);
save_figure('task2_roots.png')

contourf(Kp_x, Kd_x, z, 'ShowText', 'on');
colorbar;
xlabel('Kp_x');
ylabel('Kd_x');
title('Real Part of Roots of the Closed-Loop System', 'Interpreter', 'latex');
save_figure('task2_contour.png')

% Use a global PD controller with the gains determined from the plot
Kp_x = -79;
Kd_x = -46;

pd = @(x) Kp_x * x(1) + Kd_x * x(2) + Kp_theta * x(3) + Kd_theta * x(4); % PD control law

odefun = @(t, x) f(x, pd(r(t) - x), d(t));

[t, x] = ode45(odefun, tspan, x0);

pos = x(:, 1);
theta = x(:, 3) .* 180 / pi;
theta_dot = x(:, 4);

control = zeros(length(t), 1);

for i = 1:length(t)
    error = r(t(i)) - x(i, :)';
    control(i) = pd(error);
end

error = r(t) - x(:, :)';
max_pos = max(abs(pos));
max_theta = max(abs(theta));
max_control = max(abs(control));

figure;
plot(t, theta, 'DisplayName', 'Pendulum Angle', 'LineWidth', 1.5);
ylim([-max_theta max_theta]);
ylabel('Pendulum Angle [°]');
hold on;
yyaxis right;
plot(t, pos, 'DisplayName', 'Cart Position', 'LineWidth', 1.5);
ylim([-max_pos max_pos]);
ylabel('Cart Position [m]');
xlabel('Time [s]');
grid on;
legend('$\theta(t)$', '$x(t)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Pendulum Angle and Cart Position - Global PD Control', 'Interpreter', 'latex');

save_figure('task2_control_global.png')
