clc; clear;

addpath(genpath('.'));
addpath(genpath('../lib'));

params.density = 1000; % [kg/m^3]
params.h = 1; % [m]
params.d = 1; % [m]
params.g = 9.81; % [m/s^2]

damping = 0.003;

%% Task 1: According to an equivalent mechanical model
% of the sloshing dynamics of the cylindrical tank,
% show the formulas and the corresponding
% numerical values of the parameters (masses
% mn, lengths ln, fixed mass m0, natural
% frequencies 𝜔n(rad/s) and fn (Hz)) of the
% multi-pendulum model encompassing the
% first 10 slosh modes (n = 1,…,10).

pendulums = sloshing_pendulums(params, 10);

%% Task 2a: Derive the set of linearized EOM and the corresponding
% state-space formulation for the multi-pendulum model
% of order 10 (n = 10) having as input a prescribed lateral
% acceleration of the tank and as output the net force F_x
% exerted on the tank along x-direction.

sys_undamped = sloshing_undamped(pendulums, params);

figure;
impulse(sys_undamped, 10);
grid on;
title('Impulse response of the undamped system, n = 10');
xlabel('Time [s]');
ylabel('$F_x$ [N]', 'Interpreter', 'latex');

save_figure('task1_undamped.png')

%% Task 2b: Include in the previous EOM and related state-space
% model a modal damping ratio 𝛾n for each component
% of the multi-pendulum system.

sys_damped = sloshing_damped(pendulums, params, damping);

figure;
impulse(sys_damped, 10);
grid on;
title('Impulse response of the damped system, n = 10');
xlabel('Time [s]');
ylabel('$F_x$ [N]', 'Interpreter', 'latex');

save_figure('task1_damped.png')

%% Task 3: Using a reduced-order damped mechanical model
% including only the first fundamental slosh mode (first
% order model, n = 1), compute the time history from 0 to
% 10 s of the net force Fx due to an impulse acceleration
% of the tank using:
% 1. the built-in MATLAB® function “impulse”;
% 2. the analytical solution;
% 3. the analytical solution in modal form;
% 4. a numerical integration technique.

t_i = 0;
t_f = 10;
t_step = 0.01;
t_intervals = t_i:t_step:t_f;

pendulums = sloshing_pendulums(params, 1);
sys_damped = sloshing_damped(pendulums, params, damping);

% 1: built-in function "impulse"
[y_impulse, t_impulse] = impulse(sys_damped, t_intervals);

% 2: analytical solution
f_analitycal = @(t) sys_damped.C * expm(sys_damped.A * t) * sys_damped.B;
t_analitycal = 0:0.01:t_f;
y_analitycal = arrayfun(f_analitycal, t_analitycal);

% 3: analytical solution in modal form
[eig_vectors_A, eig_values_A] = eig(sys_damped.A);
V = eig_vectors_A;
lambda = eig_values_A;
V_inv = V \ eye(size(V));

sys_modal_A = V_inv * sys_damped.A * V;
sys_modal_B = V_inv * sys_damped.B;
sys_modal_C = sys_damped.C * V;
sys_modal_D = sys_damped.D;
sys_modal = ss(sys_modal_A, sys_modal_B, sys_modal_C, sys_modal_D);

f_modal = @(t) sys_modal.C * diag(exp(diag(lambda * t))) * sys_modal.B;

y_modal = arrayfun(f_modal, t_intervals);

% 4: numerical integration technique

% The impulse force is applied through the initial condition
x0 = [0; -1 / pendulums.L(1)];

% Comment how the impulse contributes to the initial condition, through conservation of momentum
odefun = @(t, x) sys_damped.A * x;

[t_numerical, x_numerical] = ode45(odefun, [0, t_f], x0);
y_numerical = sys_damped.C * x_numerical';

figure;
plot(t_intervals, y_impulse, 'DisplayName', 'Impulse');
hold on;
plot(t_intervals, y_analitycal, 'DisplayName', 'Analitycal');
plot(t_intervals, y_modal, 'DisplayName', 'Modal');
plot(t_numerical, y_numerical, 'DisplayName', 'Numerical');
grid on;
title('Impulse response of the damped system, n = 1');
xlabel('Time [s]');
ylabel('$F_x$ [N]', 'Interpreter', 'latex');
hold off;
legend;

save_figure('task3_impulse.png')

% Compare using ode23, ode45 and ode89
t_numerical_45 = t_numerical;
y_numerical_45 = y_numerical;
[t_numerical_23, x_numerical_23] = ode23(odefun, [0, t_f], x0);
y_numerical_23 = sys_damped.C * x_numerical_23';
[t_numerical_89, x_numerical_89] = ode89(odefun, [0, t_f], x0);
y_numerical_89 = sys_damped.C * x_numerical_89';

ode23_errors = arrayfun(@(i) abs(y_numerical_23(i) - f_modal(t_numerical_23(i))), 1:length(t_numerical_23));
ode45_errors = arrayfun(@(i) abs(y_numerical_45(i) - f_modal(t_numerical_45(i))), 1:length(t_numerical_45));
ode89_errors = arrayfun(@(i) abs(y_numerical_89(i) - f_modal(t_numerical_89(i))), 1:length(t_numerical_89));

figure;
plot(t_numerical_23, ode23_errors, 'DisplayName', 'ode23', 'LineWidth', 1.5);
hold on;
plot(t_numerical_45, ode45_errors, 'DisplayName', 'ode45', 'LineWidth', 1.5);
plot(t_numerical_89, ode89_errors, 'DisplayName', 'ode89', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Error [N]');
title('Error of the numerical integration methods');
legend('Location', 'best');
hold off;

save_figure('task3_errors.png')

%% Task 4: Find the required order n of the damped multi-
% pendulum model to achieve a reasonable desired
% accuracy on the evaluation of the net force Fx due to
% an impulse acceleration of the tank.

n_options = 1:10;
n_reference = 100;
y = zeros(length(n_options), length(t_intervals));

for n = n_options
    pendulums = sloshing_pendulums(params, n);
    sys_damped = sloshing_damped(pendulums, params, damping);

    % Use modal form to compute the impulse response
    [eig_vectors_A, eig_values_A] = eig(sys_damped.A);
    V = eig_vectors_A;
    lambda = eig_values_A;
    V_inv = V \ eye(size(V));

    y(n, :) = arrayfun(@(t) real(sys_damped.C * V * diag(exp(diag(lambda * t))) * V_inv * sys_damped.B), t_intervals);
end

pendulums = sloshing_pendulums(params, n_reference);
sys_damped = sloshing_damped(pendulums, params, damping);

y_reference = impulse(sys_damped, t_intervals)';

max_y = max(abs(y_reference));

avg_errors = zeros(length(n_options), 1);
max_errors = zeros(length(n_options), 1);

for i = 1:length(n_options)
    avg_errors(i) = mean(abs(y(i, :) - y_reference)) / max_y;
    max_errors(i) = max(abs(y(i, :) - y_reference)) / max_y;
end

figure;
plot(n_options(1:end), 100 .* avg_errors, 'DisplayName', 'Average error');
hold on;
plot(n_options(1:end), 100 .* max_errors, 'DisplayName', 'Max error');

% Trace a line at 1% error
yline(1, 'r--', '1% error', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LabelOrientation', 'horizontal', 'DisplayName', '1% error');

grid on;
xlabel('Order n of the model');
ylabel('Error [%]');
title('Relative error with respect to the maximum value of the reference model');
legend('Location', 'best');
hold off;

save_figure('task4_error.png')

figure;
plot(t_intervals, y(1, :), 'DisplayName', 'n = 1');
hold on;

for i = 2:length(n_options)
    plot(t_intervals, y(i, :), 'DisplayName', ['n = ', num2str(i)]);
end

grid on;
xlabel('Time [s]');
ylabel('Impulse response [N]');
title('Impulse response of the damped system for different orders n');
legend('Location', 'best');
hold off;

save_figure('task4_impulse.png')

%% Task 5: Using the damped multi-pendulum model of order n
% obtained from the analysis of Task 4, compute the
% time history from 0 to 100 s of the net force Fx due to a
% step acceleration of the tank using:
% 1. the built-in MATLAB® function “step”;
% 2. the analytical solution;
% 3. a numerical integration technique.

t_f = 100;

n = 5; % Order of the model
pendulums = sloshing_pendulums(params, n);
sys_damped = sloshing_damped(pendulums, params, damping);

% 1: built-in function "step"
[y_step, t_step] = step(sys_damped, t_f);

% 2: analytical solution

u = @(t) 1;

t_step_analitycal = 0:0.01:t_f;
C_inv_A = sys_damped.C / sys_damped.A;

f_analitycal = @(t) C_inv_A * (expm(sys_damped.A * t) - eye(size(sys_damped.A))) * sys_damped.B + sys_damped.D * u(t);
y_step_analytical = arrayfun(f_analitycal, t_step_analitycal);

% 3: numerical integration technique
x0 = zeros(2 * n, 1); % Initial condition for the state vector

odefun = @(t, x) sys_damped.A * x + sys_damped.B;
[t_numerical_45, x_numerical_45] = ode45(odefun, [0, t_f], x0);
y_numerical_45 = sys_damped.C * x_numerical_45' + sys_damped.D * u(t_numerical_45);

[t_numerical_23, x_numerical_23] = ode23(odefun, [0, t_f], x0);
y_numerical_23 = sys_damped.C * x_numerical_23' + sys_damped.D * u(t_numerical_23);

% Extra: add response with "frozen liquid"
t_frozen = 0:0.1:t_f;
total_mass = -sum(pendulums.m) - pendulums.m0;
y_frozen = arrayfun(@(t) total_mass * u(t), t_frozen);

figure;
plot(t_step, y_step, 'DisplayName', 'Step', 'LineWidth', 1.5);
hold on;
plot(t_step_analitycal, y_step_analytical, 'DisplayName', 'Analitycal', 'LineWidth', 1.5);
plot(t_numerical_23, y_numerical_23, 'DisplayName', 'Numerical: ode23', 'LineWidth', 1.5);
plot(t_numerical_45, y_numerical_45, 'DisplayName', 'Numerical: ode45', 'LineWidth', 1.5);
plot(t_frozen, y_frozen, 'DisplayName', 'Frozen liquid', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Step response [N]');
title('Step response of the damped system for n = 5');
legend('Location', 'best');
hold off;

save_figure('task5_step.png')

%% Step 6: Derive the frequency response function (FRF)
% corresponding to a damped n–order multi-
% pendulum model having as input the lateral
% displacement of the tank and as output the net
% force Fx exerted on the tank along x-direction.
% Show and comment every step of the derivation.

% Compare this response with that obtained by
% considering a “frozen liquid” (no sloshing).

slosh_amplitude = @(omega) abs((sys_damped.C / (1i * omega * eye(size(sys_damped.A)) - sys_damped.A) * sys_damped.B + sys_damped.D) .* (-omega .^ 2));
slosh_phase = @(omega) angle((sys_damped.C / (1i * omega * eye(size(sys_damped.A)) - sys_damped.A) * sys_damped.B + sys_damped.D) .* (-omega .^ 2));

frozen_amplitude = @(omega) abs((-sum(pendulums.m) - pendulums.m0) .* (-omega .^ 2));
frozen_phase = @(omega) angle((-sum(pendulums.m) - pendulums.m0) .* (-omega .^ 2));

omega = logspace(-1, 2, 1000); % Frequency range for the Bode plot
slosh_amplitude_values = arrayfun(slosh_amplitude, omega);
slosh_phase_values = arrayfun(slosh_phase, omega);

frozen_amplitude_values = arrayfun(frozen_amplitude, omega);
frozen_phase_values = arrayfun(frozen_phase, omega);

figure;
semilogx(omega, 20 * log10(slosh_amplitude_values), 'DisplayName', 'Amplitude', 'LineWidth', 1.5);
hold on;
semilogx(omega, 20 * log10(frozen_amplitude_values), 'DisplayName', 'Frozen liquid', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('Frequency [rad/s]');
ylabel('Magnitude [dB]');
title('Frequency response of the damped system, n = 5');
legend('Location', 'best');

save_figure('task6_amplitude.png')

figure;
semilogx(omega, slosh_phase_values * 180 / pi, 'DisplayName', 'Phase', 'LineWidth', 1.5);
hold on;
semilogx(omega, frozen_phase_values * 180 / pi, 'DisplayName', 'Frozen liquid', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('Frequency [rad/s]');
ylabel('Phase [degrees]');
title('Frequency response of the damped system, n = 5');
legend('Location', 'best');

save_figure('task6_phase.png')

%% Task 7: Step Response using SIMULINK

n = 5; % Order of the model
pendulums = sloshing_pendulums(params, n);
sys_damped = sloshing_damped(pendulums, params, damping);

x0 = zeros(2 * n, 1); % Initial condition for the state vector
u = @(t) 1;

C_inv_A = sys_damped.C / sys_damped.A;
f_analitycal = @(t) C_inv_A * (expm(sys_damped.A * t) - eye(size(sys_damped.A))) * sys_damped.B + sys_damped.D * u(t);

sim_step_ode23 = sim('step_response_ode23.slx');
sim_step_ode45 = sim('step_response_ode45.slx');

u_ode23 = sim_step_ode23.u; % Forcing term
u_ode45 = sim_step_ode45.u;

F_ode23 = sim_step_ode23.F; % Timeseries
F_ode45 = sim_step_ode45.F;

t_ode23 = F_ode23.Time; % Time vector
t_ode45 = F_ode45.Time;

y_ode23 = F_ode23.Data; % Output of the system
y_ode45 = F_ode45.Data;

dt_ode23 = diff(t_ode23); % Time intervals
dt_ode45 = diff(t_ode45);

dt_ode23_mov50 = movmean(dt_ode23, 50); % 50-point moving average
dt_ode45_mov50 = movmean(dt_ode45, 50);

ode23_error = arrayfun(@(i) abs(y_ode23(i) - f_analitycal(t_ode23(i))), 1:length(t_ode23));
ode45_error = arrayfun(@(i) abs(y_ode45(i) - f_analitycal(t_ode45(i))), 1:length(t_ode45));

ode23_error_mov50 = movmean(ode23_error, 50); % 50-point moving average
ode45_error_mov50 = movmean(ode45_error, 50); %

t_f = F_ode23.Time(end); % Final time for the step response
[y_step, t_step] = step(sys_damped, t_f); % Simulink step response

std_color = '#292f36';
ode23_color = '#20639b';
ode23_color_accent = '#3caea3';
ode45_color = '#fc915e';
ode45_color_accent = '#d53e4f';

% ode23 solution compared to the analytical solution
figure;
subplot(1, 2, 1);
plot(t_ode23, y_ode23, 'DisplayName', 'Numerical: ode23', 'LineWidth', 1.5, 'Color', std_color);
hold on;
plot(t_step', y_step, 'DisplayName', 'Step', 'LineWidth', 1.5, 'Color', ode23_color);
grid on;
xlabel('Time [s]');
ylabel('Step response [N]');
title('Step response - \texttt{ode23}', 'Interpreter', 'latex');
legend('Location', 'best');
hold off;

% ode45 solution compared to the analytical solution
subplot(1, 2, 2);
plot(t_ode45, y_ode45, 'DisplayName', 'Numerical: ode45', 'LineWidth', 1.5, 'Color', std_color);
hold on;
plot(t_step, y_step, 'DisplayName', 'Step', 'LineWidth', 1.5, 'Color', ode45_color);
grid on;
xlabel('Time [s]');
ylabel('Step response [N]');
title('Step response - \texttt{ode45}', 'Interpreter', 'latex');
legend('Location', 'best');
hold off;

save_figure('task7_step_composite.png', keep_title = true, aspect_ratio_multiplier = 1.5)

% ode23 time intervals
figure;
subplot(1, 2, 1);
plot(dt_ode23, 'DisplayName', 'Time intervals', 'LineWidth', 1.5, 'Color', ode23_color);
hold on;
plot(dt_ode23_mov50, 'DisplayName', '50-point moving average', 'LineWidth', 1.5, 'Color', ode23_color_accent);
grid on;
xlabel('Index');
ylabel('Time interval [s]');
title('Time intervals - \texttt{ode23}', 'Interpreter', 'latex');
legend('Location', 'best');
hold off;

% ode45 time intervals
subplot(1, 2, 2);
plot(dt_ode45, 'DisplayName', 'Time intervals', 'LineWidth', 1.5, 'Color', ode45_color);
hold on;
plot(dt_ode45_mov50, 'DisplayName', '50-point moving average', 'LineWidth', 1.5, 'Color', ode45_color_accent);
grid on;
xlabel('Index');
ylabel('Time interval [s]');
title('Time intervals - \texttt{ode45}', 'Interpreter', 'latex');
legend('Location', 'best');
hold off;

save_figure('task7_time_intervals_composite.png', keep_title = true, aspect_ratio_multiplier = 1.5)

% ode23 error compared to the analytical solution
figure;
subplot(1, 2, 1);
plot(t_ode23, ode23_error, 'DisplayName', 'Error', 'LineWidth', 1.5, 'Color', ode23_color);
hold on;
plot(t_ode23, ode23_error_mov50, 'DisplayName', '50-point moving average', 'LineWidth', 1.5, 'Color', ode23_color_accent);
grid on;
xlabel('Time [s]');
ylabel('Error [N]');
title('Error - \texttt{ode23}', 'Interpreter', 'latex');
legend('Location', 'best');
hold off;

% ode45 error compared to the analytical
subplot(1, 2, 2);
plot(t_ode45, ode45_error, 'DisplayName', 'Error', 'LineWidth', 1.5, 'Color', ode45_color);
hold on;
plot(t_ode45, ode45_error_mov50, 'DisplayName', '50-point moving average', 'LineWidth', 1.5, 'Color', ode45_color_accent);
grid on;
xlabel('Time [s]');
ylabel('Error [N]');
title('Error - \texttt{ode45}', 'Interpreter', 'latex');
legend('Location', 'best');
hold off;

save_figure('task7_error_composite.png', keep_title = true, aspect_ratio_multiplier = 1.5)
