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

%% Task 2b: Include in the previous EOM and related state-space
% model a modal damping ratio 𝛾n for each component
% of the multi-pendulum system.

sys_damped = sloshing_damped(pendulums, params, damping);

figure;
impulse(sys_damped, 10);
grid on;

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

% s * X = A * X + B * U
% (s * I - A) * X = B * U
% X = (s * I - A)^-1 * B * U
% Y = C * X + D * U = C * (s * I - A)^-1 * B * U + D * U

syms s t;
U = 1;
Y = simplify(sys_damped.C / (s * eye(size(sys_damped.A)) - sys_damped.A) * sys_damped.B * U + sys_damped.D * U);
y = matlabFunction(ilaplace(Y, s, t));

t_analitycal = 0:0.01:t_f;
y_analitycal = y(t_analitycal);

% 3: TODO: analytical solution in modal form

% 4: numerical integration technique

% The impulse force is applied through the initial condition
x0 = zeros(size(sys_damped.A, 1), 1);
% TODO: find the initial condition for the impulse force
x0(2) = 1;

odefun = @(t, x) sys_damped.A * x;
[t_numerical, x_numerical] = ode45(odefun, [0, t_f], x0);
y_numerical = sys_damped.C * x_numerical';

figure;
plot(t_impulse, y_impulse, 'DisplayName', 'impulse');
hold on;
plot(t_analitycal, y_analitycal, 'DisplayName', 'analitycal');
plot(t_numerical, y_numerical, 'DisplayName', 'numerical');
grid on;
legend;
