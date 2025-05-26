clc; clear;

addpath(genpath('.'));
addpath(genpath('../lib'));

% Aircraft parameters
params.m = 7484.4; % mass [Kg]
params.I_yy = 84309; % pitching moment of inertia [Kg*m^2]
params.S = 32.8; % wing area [m^2]
params.c = 2.29; % mean aerodynamic chord [m]
params.a_T = 0; % thrustline angle [rad]
params.z_T = 0.378; % thrustline vertical distance [m]
params.g = 9.80665; % gravity [m/s^2]

% Get the aerodynamic model
aer = aerodynamic_model(params);

%% Task 1  – Trim problem
% Starting from the model of longitudinal dynamics in the
% body axes frame ( (u, w, q, theta, h) - model, see part 1),
% derive the equations governing the trim problem and
% determine numerically the trim state (equilibrium flight
% condition) corresponding to a steady-state straight-and-
% level flight at sea level with airspeed V = 120 knots.

V_trim = 120 * 1.852/3.6; % airspeed [m/s]
h_trim = 0; % altitude [m]

f = longitudinal_model(params, aer);

% x_trim is the solution of the trim problem
% 1: angle of attack (alpha)
% 2: elevator deflection (delta)
% 3: thrust (T)

trim_solution = fsolve(@(x) trim_eq(x, V_trim, h_trim, f), zeros(3, 1));

%% Task 2 – Trim verification
% Verify the trim state found in Task 1 by developing a MATLAB® code
% to simulate the response of the aircraft in the time interval [0-5]s
% when x(0) = x_trim; and the thrust and the elevator are held at the calculated trim values.
% Plot the time response in terms of forward speed, heave velocity,
% pitch rate, pitch attitude, altitude and angle of attack.

T_trim = trim_solution(1);
alpha_trim = trim_solution(2);
delta_trim = trim_solution(3);
x_trim = [
          V_trim * cos(alpha_trim); % u [m/s]
          V_trim * sin(alpha_trim); % w [m/s]
          0; % q [rad/s]
          alpha_trim; % theta [rad]
          h_trim % h [m]
          ];
u_trim = [delta_trim; T_trim];

u = @(t) u_trim;
tspan = [0 20];

odefun = @(t, x) f(x, u(t));

[t, x] = ode45(odefun, tspan, x_trim);
u = x(:, 1);
w = x(:, 2);
q_deg = x(:, 3) * 180 / pi;
theta_deg = x(:, 4) * 180 / pi;
h = x(:, 5);

alpha = atan(w ./ u);
alpha_deg = alpha * 180 / pi;

figure;
subplot(3, 2, 1);
plot(t, u);
grid on;
xlabel('Time [s]');
ylabel('$u$ [m/s]', 'Interpreter', 'latex');
title('Forward speed', 'Interpreter', 'latex');

subplot(3, 2, 2);
plot(t, w);
grid on;
xlabel('Time [s]');
ylabel('$w$ [m/s]', 'Interpreter', 'latex');
title('Heave velocity', 'Interpreter', 'latex');

subplot(3, 2, 3);
plot(t, q_deg);
grid on;
xlabel('Time [s]');
ylabel('q [deg/s]', 'Interpreter', 'latex');
title('Pitch rate', 'Interpreter', 'latex');

subplot(3, 2, 4);
plot(t, theta_deg);
grid on;
xlabel('Time [s]');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
title('Pitch attitude', 'Interpreter', 'latex');

subplot(3, 2, 5);
plot(t, h);
grid on;
xlabel('Time [s]');
ylabel('$h$ [m]', 'Interpreter', 'latex');
title('Altitude', 'Interpreter', 'latex');

subplot(3, 2, 6);
plot(t, alpha_deg);
grid on;
xlabel('Time [s]');
ylabel('$\alpha$ [deg]', 'Interpreter', 'latex');
title('Angle of attack', 'Interpreter', 'latex');

sgtitle('Longitudinal dynamics', 'Interpreter', 'latex');

save_figure('task2_trim.png', keep_title = true);

figure;
subplot(3, 2, 1);
plot(t, u);
ylim([60, 63]);
grid on;
xlabel('Time [s]');
ylabel('$u$ [m/s]', 'Interpreter', 'latex');
title('Forward speed', 'Interpreter', 'latex');

subplot(3, 2, 2);
plot(t, w);
ylim([0.85, 0.95]);
grid on;
xlabel('Time [s]');
ylabel('$w$ [m/s]', 'Interpreter', 'latex');
title('Heave velocity', 'Interpreter', 'latex');

subplot(3, 2, 3);
plot(t, q_deg);
ylim([-0.1, 0.1]);
grid on;
xlabel('Time [s]');
ylabel('q [deg/s]', 'Interpreter', 'latex');
title('Pitch rate', 'Interpreter', 'latex');

subplot(3, 2, 4);
plot(t, theta_deg);
ylim([0.8, 0.9]);
grid on;
xlabel('Time [s]');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
title('Pitch attitude', 'Interpreter', 'latex');

subplot(3, 2, 5);
plot(t, h);
ylim([-0.1, 0.1]);
grid on;
xlabel('Time [s]');
ylabel('$h$ [m]', 'Interpreter', 'latex');
title('Altitude', 'Interpreter', 'latex');

subplot(3, 2, 6);
plot(t, alpha_deg);
ylim([0.8, 0.9]);
grid on;
xlabel('Time [s]');
ylabel('$\alpha$ [deg]', 'Interpreter', 'latex');
title('Angle of attack', 'Interpreter', 'latex');

sgtitle('Longitudinal dynamics', 'Interpreter', 'latex');

save_figure('task2_trim_rescaled.png', keep_title = true);

%% Task 3 - Nonlinear response to elevator pulse
% Starting from the initial condition corresponding to x(0) = x_trim,
% develop a MATLAB® code to simulate the nonlinear response of the
% aircraft in the time interval [0-100]s, when the thrust is held at the
% calculated trim value and the elevator is subjected to a one-degree
% pulse deflection from the trim value in the interval [5-10]s:
%
% delta(t) = delta_trim + Delta delta, Delta delta = -1 deg.
%
% A. Plot the time response in terms of forward speed, heave velocity,
% pitch rate, pitch attitude, altitude and angle of attack.
%
% B. Provide comments on the evolution of the variables.

delta = @(t) delta_trim + (t >= 5 & t <= 10) * (-1 * pi / 180); % elevator deflection [rad]
input = @(t) [delta(t); T_trim];

tspan = [0 100];
odefun = @(t, x) f(x, input(t));
[t, x] = ode45(odefun, tspan, x_trim);

u = x(:, 1);
w = x(:, 2);
q_deg = x(:, 3) * 180 / pi;
theta_deg = x(:, 4) * 180 / pi;
h = x(:, 5);
alpha = atan(w ./ u);
alpha_deg = alpha * 180 / pi;

figure;
subplot(3, 2, 1);
plot(t, u);
grid on;
xlabel('Time [s]');
ylabel('$u$ [m/s]', 'Interpreter', 'latex');
title('Forward speed', 'Interpreter', 'latex');

subplot(3, 2, 2);
plot(t, w);
grid on;
xlabel('Time [s]');
ylabel('$w$ [m/s]', 'Interpreter', 'latex');
title('Heave velocity', 'Interpreter', 'latex');

subplot(3, 2, 3);
plot(t, q_deg);
grid on;
xlabel('Time [s]');
ylabel('$q$ [deg/s]', 'Interpreter', 'latex');
title('Pitch rate', 'Interpreter', 'latex');

subplot(3, 2, 4);
plot(t, theta_deg);
grid on;
xlabel('Time [s]');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
title('Pitch attitude', 'Interpreter', 'latex');

subplot(3, 2, 5);
plot(t, h);
grid on;
xlabel('Time [s]');
ylabel('$h$ [m]', 'Interpreter', 'latex');
title('Altitude', 'Interpreter', 'latex');

subplot(3, 2, 6);
plot(t, alpha_deg);
grid on;
xlabel('Time [s]');
ylabel('$\alpha$ [deg]', 'Interpreter', 'latex');
title('Angle of attack', 'Interpreter', 'latex');

sgtitle('Longitudinal dynamics - Longitudinal pulse response', 'Interpreter', 'latex');

save_figure('task3_nonlinear_response.png', keep_title = true);

%% Task 4 – Nonlinear response to elevator pulse (Simulink)
% Repeat Task 3 using a Simulink model.

sim_nl = sim("task4_simulink.slx");

t = sim_nl.tout;

u = sim_nl.x.signals.values(:, 1, :);
w = sim_nl.x.signals.values(:, 2, :);
q_deg = sim_nl.x.signals.values(:, 3, :) * 180 / pi;
theta_deg = sim_nl.x.signals.values(:, 4, :) * 180 / pi;
h = sim_nl.x.signals.values(:, 5, :);
alpha = atan(w ./ u);
alpha_deg = alpha * 180 / pi;

figure;
subplot(3, 2, 1);
plot(t, u);
grid on;
xlabel('Time [s]');
ylabel('$u$ [m/s]', 'Interpreter', 'latex');
title('Forward speed', 'Interpreter', 'latex');

subplot(3, 2, 2);
plot(t, w);
grid on;
xlabel('Time [s]');
ylabel('$w$ [m/s]', 'Interpreter', 'latex');
title('Heave velocity', 'Interpreter', 'latex');

subplot(3, 2, 3);
plot(t, q_deg);
grid on;
xlabel('Time [s]');
ylabel('$q$ [deg/s]', 'Interpreter', 'latex');
title('Pitch rate', 'Interpreter', 'latex');

subplot(3, 2, 4);
plot(t, theta_deg);
grid on;
xlabel('Time [s]');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
title('Pitch attitude', 'Interpreter', 'latex');

subplot(3, 2, 5);
plot(t, h);
grid on;
xlabel('Time [s]');
ylabel('$h$ [m]', 'Interpreter', 'latex');
title('Altitude', 'Interpreter', 'latex');

subplot(3, 2, 6);
plot(t, alpha_deg);
grid on;
xlabel('Time [s]');
ylabel('$\alpha$ [deg]', 'Interpreter', 'latex');
title('Angle of attack', 'Interpreter', 'latex');

sgtitle('Longitudinal dynamics - Longitudinal pulse response - Simulink', 'Interpreter', 'latex');

save_figure('task4_simulink_response.png', keep_title = true);

%% Task 5 – Linearized model – EOM
% stability.X_u = -0.057076461;
% stability.X_w = 0.125051144;
% stability.X_q = 0;
% stability.X_w_dot = 0;
% stability.Z_u = -0.30507998;
% stability.Z_w = -0.86334716;
% stability.Z_q = 0;
% stability.Z_w_dot = 0;
% stability.M_u = -0.0014736;
% stability.M_w = -0.036668742;
% stability.M_q = -0.5442436;
% stability.M_w_dot = 0;
% stability.X_delta = 0;
% stability.Z_delta = -7.38503134;
% stability.M_delta = -3.90965145;

% stability.X_T = 0;
% stability.Z_T = 0;
% stability.M_T = 0;

stability = longitudinal_derivatives(params, aer, x_trim, u_trim);

[A, B, C, D] = longitudinal_linear_model(params, stability, x_trim);

%% Task 6.1 – Analysis of the linearized system
% By using the numerical values of stability and
% control derivatives shown in the table alongside,
% report the natural frequency and damping ratio
% associated with the four eigenvalues of the
% linearized system.
% Provide comments on the effect of their values on
% the longitudinal dynamic characteristics.

eig_A = eig(A);

min_real = min(real(eig_A));
max_real = max(real(eig_A));
min_imag = min(imag(eig_A));
max_imag = max(imag(eig_A));

figure;
plot(real(eig_A), imag(eig_A), "x", 'MarkerSize', 15);
xlim([min_real - 0.3, max_real + 0.3]);
ylim([min_imag - 0.3, max_imag + 0.3]);
grid on;
title("Eigenvalues of State Matrix A", Interpreter = 'latex');

save_figure('task6_eigenvalues.png', keep_title = true);

w_n = abs(eig_A);
f_n = 2 * pi * w_n;

xi_n = cos(atan(imag(eig_A) ./ real(eig_A)));

fprintf("Eigenvalues of A:\n");
disp(eig_A);
fprintf("Natural frequencies (rad/s):\n");
disp(w_n);
fprintf("Frequencies (Hz):\n");
disp(f_n);
fprintf("Damping ratios:\n");
disp(xi_n);

%% Task 6.2 – Analysis of the linearized system
% Considering the elevator deflection as the control input,
% plot the Bode diagram (magnitude and phase) and comment
% the resulting behaviour of the following transfer functions:

f_amplitude = @(omega) abs(C / (1i * omega * eye(size(A)) - A) * B + D);
f_phase = @(omega) angle(C / (1i * omega * eye(size(A)) - A) * B + D);

omega = logspace(-2, 3, 1000);

G_u_amplitude = arrayfun(@(omega) at(f_amplitude(omega), 1, 1), omega);
G_w_amplitude = arrayfun(@(omega) at(f_amplitude(omega), 2, 1), omega);
G_q_amplitude = arrayfun(@(omega) at(f_amplitude(omega), 3, 1), omega);
G_theta_amplitude = arrayfun(@(omega) at(f_amplitude(omega), 4, 1), omega);

phase_clock = @(x) (x > 0) .* (x - 2 .* pi) + (x <= 0) .* x;

G_u_phase = phase_clock(arrayfun(@(omega) at(f_phase(omega), 1, 1), omega)) .* 180 / pi;
G_w_phase = phase_clock(arrayfun(@(omega) at(f_phase(omega), 2, 1), omega)) .* 180 / pi;
G_q_phase = phase_clock(arrayfun(@(omega) at(f_phase(omega), 3, 1), omega)) .* 180 / pi;
G_theta_phase = phase_clock(arrayfun(@(omega) at(f_phase(omega), 4, 1), omega)) .* 180 / pi;

db = @(x) 20 * log10(x);

[G_u_a, G_u_b] = ss2tf(A, B(:, 1), C(1, :), D(1, 1));
[G_w_a, G_w_b] = ss2tf(A, B(:, 1), C(2, :), D(2, 1));
[G_q_a, G_q_b] = ss2tf(A, B(:, 1), C(3, :), D(3, 1));
[G_theta_a, G_theta_b] = ss2tf(A, B(:, 1), C(4, :), D(4, 1));

G_u = tf(G_u_a, G_u_b);
G_w = tf(G_w_a, G_w_b);
G_q = tf(G_q_a, G_q_b);
G_theta = tf(G_theta_a, G_theta_b);

[G_u_amplitude_asymp, G_u_phase_asymp] = asymp_bode(G_u, omega);
[G_w_amplitude_asymp, G_w_phase_asymp] = asymp_bode(G_w, omega);
[G_q_amplitude_asymp, G_q_phase_asymp] = asymp_bode(G_q, omega);
[G_theta_amplitude_asymp, G_theta_phase_asymp] = asymp_bode(G_theta, omega);

figure;
subplot(2, 1, 1);
semilogx(omega, db(G_u_amplitude), 'LineWidth', 1.5);
hold on;
semilogx(omega, db(G_u_amplitude_asymp), 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
legend('$|G_{u\delta_e}(j\omega)|$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_{u\delta_e}$', 'Interpreter', 'latex');

subplot(2, 1, 2);
semilogx(omega, G_u_phase, 'LineWidth', 1.5);
hold on;
semilogx(omega, G_u_phase_asymp, 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [°]');
legend('$\angle G_{u\delta_e}(j\omega)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_{u\delta_e}(j\omega)$', 'Interpreter', 'latex');

save_figure('task6_bode_G_u.png', keep_title = true);

figure;
subplot(2, 1, 1);
semilogx(omega, db(G_w_amplitude), 'LineWidth', 1.5);
hold on;
semilogx(omega, db(G_w_amplitude_asymp), 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
legend('$|G_{w\delta_e}(j\omega)|$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_{w\delta_e}$', 'Interpreter', 'latex');

subplot(2, 1, 2);
semilogx(omega, G_w_phase, 'LineWidth', 1.5);
hold on;
semilogx(omega, G_w_phase_asymp, 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [°]');
legend('$\angle G_{w\delta_e}(j\omega)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_{w\delta_e}(j\omega)$', 'Interpreter', 'latex');

save_figure('task6_bode_G_w.png', keep_title = true);

figure;
subplot(2, 1, 1);
semilogx(omega, db(G_q_amplitude), 'LineWidth', 1.5);
hold on;
semilogx(omega, db(G_q_amplitude_asymp), 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
legend('$|G_{q\delta_e}(j\omega)|$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_{q\delta_e}$', 'Interpreter', 'latex');

subplot(2, 1, 2);
semilogx(omega, G_q_phase, 'LineWidth', 1.5);
hold on;
semilogx(omega, G_q_phase_asymp, 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [°]');
legend('$\angle G_{q\delta_e}(j\omega)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_{q\delta_e}(j\omega)$', 'Interpreter', 'latex');

save_figure('task6_bode_G_q.png', keep_title = true);

figure;
subplot(2, 1, 1);
semilogx(omega, db(G_theta_amplitude), 'LineWidth', 1.5);
hold on;
semilogx(omega, db(G_theta_amplitude_asymp), 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
legend('$|G_{\theta\delta_e}(j\omega)|$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_{\theta\delta_e}$', 'Interpreter', 'latex');

subplot(2, 1, 2);
semilogx(omega, G_theta_phase, 'LineWidth', 1.5);
hold on;
semilogx(omega, G_theta_phase_asymp, 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [°]');
legend('$\angle G_{\theta\delta_e}(j\omega)$', 'Interpreter', 'latex', 'Location', 'Best');
title('Bode Diagram of the Linearized System - $G_{\theta\delta_e}(j\omega)$', 'Interpreter', 'latex');

save_figure('task6_bode_G_theta.png', keep_title = true);

%% Task 7 – Linear response to elevator pulse

x0 = zeros(4, 1);
delta_delta = @(t) (t >= 5 & t <= 10) * (-1 * pi / 180); % elevator deflection [rad]
delta_input = @(t) [delta_delta(t); 0];

tspan = [0 100]; % time span
odefun = @(t, delta_x) A * delta_x + B * delta_input(t);
[t_lin, delta_x] = ode45(odefun, tspan, x0); % solve the ODE
delta_y = arrayfun(@(idx) C * delta_x(idx, :)' + D * delta_input(t_lin(idx)), 1:length(t_lin), 'UniformOutput', false);
delta_y = cell2mat(delta_y)';

u_lin = delta_y(:, 1) + x_trim(1);
w_lin = delta_y(:, 2) + x_trim(2);
q_lin = delta_y(:, 3) + x_trim(3);
q_deg_lin = q_lin * 180 / pi;
theta_lin = delta_y(:, 4) + x_trim(4);
theta_deg_lin = theta_lin * 180 / pi;
alpha_lin = delta_y(:, 6) + x_trim(4);
alpha_deg_lin = alpha_lin * 180 / pi;

delta_h_dot_lin = delta_y(:, 5);
[~, delta_h_lin] = ode45(@(t, ~) interp1(t_lin, delta_h_dot_lin, t, 'spline'), t_lin, 0);
h_lin = delta_h_lin + x_trim(5);

figure;
subplot(3, 2, 1);
plot(t, u, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, u_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$u$ [m/s]', 'Interpreter', 'latex');
title('Forward speed', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 2);
plot(t, w, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, w_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$w$ [m/s]', 'Interpreter', 'latex');
title('Heave velocity', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 3);
plot(t, q_deg, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, q_deg_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$q$ [deg/s]', 'Interpreter', 'latex');
title('Pitch rate', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 4);
plot(t, theta_deg, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, theta_deg_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
title('Pitch attitude', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 5);
plot(t, alpha_deg, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, alpha_deg_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$\alpha$ [deg]', 'Interpreter', 'latex');
title('Angle of attack', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 6);
plot(t, h, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, h_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$h$ [m]', 'Interpreter', 'latex');
title('Altitude', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

sgtitle('Longitudinal dynamics - Linearized model', 'Interpreter', 'latex');

save_figure('task7_linear_response.png', keep_title = true);

%% Task 8 - Simulink

sim_l = sim("task8_simulink.slx");

t_lin = sim_l.tout;
u_lin = sim_l.x_lin.signals.values(:, 1);
w_lin = sim_l.x_lin.signals.values(:, 2);
q_dot_lin = sim_l.x_lin.signals.values(:, 3);
theta_dot_lin = sim_l.x_lin.signals.values(:, 4);
h_lin = sim_l.h_lin;
alpha_lin = sim_l.alpha_lin;

q_deg_lin = q_dot_lin * 180 / pi;
theta_deg_lin = theta_dot_lin * 180 / pi;
alpha_deg_lin = alpha_lin * 180 / pi;

figure;
subplot(3, 2, 1);
plot(t, u, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, u_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$u$ [m/s]', 'Interpreter', 'latex');
title('Forward speed', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 2);
plot(t, w, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, w_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$w$ [m/s]', 'Interpreter', 'latex');
title('Heave velocity', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 3);
plot(t, q_deg, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, q_deg_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$q$ [deg/s]', 'Interpreter', 'latex');
title('Pitch rate', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 4);
plot(t, theta_deg, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, theta_deg_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
title('Pitch attitude', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 5);
plot(t, alpha_deg, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, alpha_deg_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$\alpha$ [deg]', 'Interpreter', 'latex');
title('Angle of attack', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

subplot(3, 2, 6);
plot(t, h, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, h_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$h$ [m]', 'Interpreter', 'latex');
title('Altitude', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

sgtitle('Longitudinal dynamics - Linearized model - Simulink', 'Interpreter', 'latex');

save_figure('task8_linear_response_simulink.png', keep_title = true);
