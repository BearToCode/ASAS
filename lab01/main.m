clc; clear;

addpath(genpath('.'));

params.density = 1000;
params.h = 1;
params.d = 1;
params.g = 9.81;

damping = 0.003;

%% Task 2a: Derive the set of linearized EOM and the corresponding
% state-space formulation for the multi-pendulum model
% of order 10 (n = 10) having as input a prescribed lateral
% acceleration of the tank and as output the net force F_x
% exerted on the tank along x-direction.

pendulums = sloshing_pendulums(params, 10);
sys_undamped = sloshing_undamped(pendulums, params);

figure;
impulse(sys_undamped, 10);
grid on;
title('Impulse response of the undamped system, n = 10');

%% Task 2b: Include in the previous EOM and related state-space
% model a modal damping ratio 𝛾n for each component
% of the multi-pendulum system.

sys_damped = sloshing_damped(pendulums, params, damping);

figure;
impulse(sys_damped, 10);
grid on;
title('Impulse response of the damped system, n = 10');

%% Task 3: Using a reduced-order damped mechanical model
% including only the first fundamental slosh mode (first
% order model, n = 1), compute the time history from 0 to
% 10 s of the net force Fx due to an impulse acceleration
% of the tank using:
% 1. the built-in MATLAB® function “impulse”;
% 2. the analytical solution;
% 3. the analytical solution in modal form;
% 4. a numerical integration technique.

t_f = 10;

pendulums = sloshing_pendulums(params, 1);
sys_damped = sloshing_damped(pendulums, params, damping);

% 1: built-in function "impulse"
[y_impulse, t_impulse] = impulse(sys_damped, t_f);

% 2: analytical solution
syms s t;
H = sys_damped.C / (s * eye(size(sys_damped.A)) - sys_damped.A) * sys_damped.B;
Y = H * 1;
f_analitycal = matlabFunction(ilaplace(Y, s, t));

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

f_modal = @(t) sys_damped.C * V * expm(lambda * t) * V_inv * sys_damped.B;

t_modal = 0:0.01:t_f;
y_modal = arrayfun(f_modal, t_modal);

% 4: numerical integration technique

% The impulse force is applied through the initial condition
x0 = [0; -1 / pendulums.L(1)];

odefun = @(t, x) sys_damped.A * x;
[t_numerical, x_numerical] = ode45(odefun, [0, t_f], x0);
y_numerical = sys_damped.C * x_numerical';

figure;
plot(t_impulse, y_impulse, 'DisplayName', 'impulse');
hold on;
plot(t_analitycal, y_analitycal, 'DisplayName', 'analitycal');
plot(t_modal, y_modal, 'DisplayName', 'modal');
plot(t_numerical, y_numerical, 'DisplayName', 'numerical');
grid on;
legend;

%% Task 4: Find the required order n of the damped multi-
% pendulum model to achieve a reasonable desired
% accuracy on the evaluation of the net force Fx due to
% an impulse acceleration of the tank.

n_options = 1:10;
t = 0:0.01:t_f;
y = zeros(length(n_options), length(t));

for n = n_options
    pendulums = sloshing_pendulums(params, n);
    sys_damped = sloshing_damped(pendulums, params, damping);

    % Use modal form to compute the impulse response
    [eig_vectors_A, eig_values_A] = eig(sys_damped.A);
    V = eig_vectors_A;
    lambda = eig_values_A;
    V_inv = V \ eye(size(V));

    y(n, :) = arrayfun(@(t) real(sys_damped.C * V * expm(lambda * t) * V_inv * sys_damped.B), t);
end

reference_y = y(end, :);
max_y = max(abs(reference_y));

avg_errors = zeros(length(n_options) - 1, 1);
max_errors = zeros(length(n_options) - 1, 1);

for i = 1:length(n_options) - 1
    avg_errors(i) = mean(abs(y(i, :) - reference_y)) / max_y;
    max_errors(i) = max(abs(y(i, :) - reference_y)) / max_y;
end

figure;
plot(n_options(1:end - 1), avg_errors, 'DisplayName', 'Average error');
hold on;
plot(n_options(1:end - 1), max_errors, 'DisplayName', 'Max error');
grid on;
xlabel('Order n of the model');
ylabel('Error [%]');
title('Error in the impulse response of the damped system for different orders n');
legend('Location', 'best');
hold off;

figure;
plot(t, y(1, :), 'DisplayName', 'n = 1');
hold on;

for i = 2:length(n_options)
    plot(t, y(i, :), 'DisplayName', ['n = ', num2str(i)]);
end

grid on;
xlabel('Time [s]');
ylabel('Impulse response [N]');
title('Impulse response of the damped system for different orders n');
legend('Location', 'best');
hold off;

%% Task 5: Using the damped multi-pendulum model of order n
% obtained from the analysis of Task 4, compute the
% time history from 0 to 100 s of the net force Fx due to a
% step acceleration of the tank using:
% 1. the built-in MATLAB® function “step”;
% 2. the analytical solution;
% 3. a numerical integration technique.

t_f = 100;

n = 4; % Order of the model
pendulums = sloshing_pendulums(params, n);
sys_damped = sloshing_damped(pendulums, params, damping);

% 1: built-in function "step"
[y_step, t_step] = step(sys_damped, t_f);

% 2: analytical solution
syms s t;
H = sys_damped.C / (s * eye(size(sys_damped.A)) - sys_damped.A) * sys_damped.B + sys_damped.D;
Y = H * 1 / s;
f_step = matlabFunction(ilaplace(Y, s, t));

t_step_analitycal = 0:0.01:t_f;
y_step_analitycal = arrayfun(f_step, t_step_analitycal);

% 3: numerical integration technique
u = @(t) 1;
x0 = zeros(2 * n, 1); % Initial condition for the state vector

odefun = @(t, x) sys_damped.A * x + sys_damped.B;
[t_numerical, x_numerical] = ode45(odefun, [0, t_f], x0);

y_numerical = sys_damped.C * x_numerical' + sys_damped.D * u(t_numerical);

figure;
plot(t_step, y_step, 'DisplayName', 'step');
hold on;
plot(t_step_analitycal, y_step_analitycal, 'DisplayName', 'analitycal');
plot(t_numerical, y_numerical, 'DisplayName', 'numerical');
grid on;
xlabel('Time [s]');
ylabel('Step response [N]');
title('Step response of the damped system for n = 4');
legend('Location', 'best');
hold off;

%% Step 6: Derive the frequency response function (FRF)
% corresponding to a damped n–order multi-
% pendulum model having as input the lateral
% displacement of the tank and as output the net
% force Fx exerted on the tank along x-direction.
% Show and comment every step of the derivation.

figure;
bode(sys_damped, {0.1, 100});
grid on;
title('Bode plot of the damped system, n = 4');

% Add the frequency response for n = 1 for reference

pendulums = sloshing_pendulums(params, 1);
sys_damped = sloshing_damped(pendulums, params, damping);

figure;
bode(sys_damped, {0.1, 100});
grid on;
title('Bode plot of the damped system, n = 1');
