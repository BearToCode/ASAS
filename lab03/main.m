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

%save_figure('task2_trim.png', keep_title = true);

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

%% Task 4 - Simulink
clc;
simnl = sim("task4_simulink.slx");
% plot(simnl.tout, simnl.delta)
figure
% sistemare nomi variabili!
for hhh = 1:5
subplot(2,3,hhh)
plotta = reshape(simnl.x.signals.values(hhh,1,:), [1, length(simnl.x.signals.values(hhh,1,:))]);
plot(simnl.tout, plotta)
end
subplot(2,3,6)
plot(simnl.tout, simnl.alpha.data)
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

[A, B] = longitudinal_linear_model(params, stability, x_trim);

C = eye(4);
D = zeros(4,2);

%% Task 6 - Rivedere
eigA = eig(A);
figure
plot(real(eigA), imag(eigA), "x")
grid on;
title("Eigenvalues of State Matrix A");

w_n = sqrt(abs(real(eigA))); % è giusto ?
f_n = 2*pi*w_n; 

xi_n = cos(atan(imag(eigA)./real(eigA))); % è giusto ?

% FUNZIONI DI TRASFERIMENTO ....................
syms s;
G = C*(s*eye(4) - A)\eye(4)*B(:,1) + D
%% Task 7 – Linear response to elevator pulse

% MANCA PLOT H ! 
x0 = zeros(4, 1);
delta_delta = @(t) (t >= 5 & t <= 10) * (-1 * pi / 180); % elevator deflection [rad]
delta_input = @(t) [delta_delta(t); 0];

tspan = [0 100]; % time span
odefun = @(t, delta_x) A * delta_x + B * delta_input(t);
[t_lin, delta_x] = ode45(odefun, tspan, x0); % solve the ODE

u_lin = delta_x(:, 1) + x_trim(1);
w_lin = delta_x(:, 2) + x_trim(2);
q_lin = delta_x(:, 3) + x_trim(3);
q_deg_lin = q_lin * 180 / pi;
theta_lin = delta_x(:, 4) + x_trim(4);
theta_deg_lin = theta_lin * 180 / pi;
alpha_lin = atan(w_lin ./ u_lin);
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

subplot(3, 2, [5 6]);
plot(t, alpha_deg, 'DisplayName', 'Nonlinear', 'LineStyle', '--');
hold on;
plot(t_lin, alpha_deg_lin, 'DisplayName', 'Linearized');
grid on;
xlabel('Time [s]');
ylabel('$\alpha$ [deg]', 'Interpreter', 'latex');
title('Angle of attack', 'Interpreter', 'latex');
legend('Interpreter', 'latex')

sgtitle('Longitudinal dynamics - Linearized model', 'Interpreter', 'latex');

save_figure('task7_linear_response.png', keep_title = true);

%% Task 8 - Simulink
% controllare
clc;
% cambiare nome a L!
L = [sin(x_trim(4)), -cos(x_trim(4)), x_trim(1)*cos(x_trim(4)) + x_trim(2)*sin(x_trim(4));...
    (-x_trim(2)/x_trim(1)^2) / (1 + (x_trim(2)/x_trim(1))^2), (1 / x_trim(1)) / (1 + (x_trim(2)/x_trim(1))^2), 0];
siml = sim("task8_simulink.slx");
% figure
% plot(siml.tout, siml.delta*180/pi)
figure
for idx = 1:4
    subplot(2,3,idx);
    plot(siml.tout, siml.y.signals.values(:, idx))
    grid on
end
subplot(2,3,5)
plot(siml.tout, siml.h)
grid on
subplot(2,3,6)
plot(siml.tout, siml.alpha)
grid on

