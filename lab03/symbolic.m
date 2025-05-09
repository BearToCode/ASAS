clc; clear;

syms S m I_yy c a_T z_T g real; % params
syms u w q theta h real; % state variables
syms w_dot real;
syms T delta real; % control inputs
syms c_L0 c_L_alpha c_L_delta c_D0 c_D_alpha c_D_alpha2 c_M0 c_M_alpha c_M_q c_M_delta real; % aerodynamics
syms rho real; % air density

V = sqrt(u ^ 2 + w ^ 2); % speed
alpha = atan(w / u); % angle of attack

c_L = c_L0 + c_L_alpha * alpha + c_L_delta * delta;
c_D = c_D0 + c_D_alpha * alpha + c_D_alpha2 * alpha ^ 2;
c_M = c_M0 + c_M_alpha * alpha + c_M_q * (c / sqrt(u ^ 2 * w ^ 2)) * q + c_M_delta * delta;

dyn_pres = 0.5 * rho * V ^ 2; % dynamic pressure
L = c_L * dyn_pres * S;
D = c_D * dyn_pres * S;

% Aerodynamics forces
X_A = L * sin(alpha) - D * cos(alpha);
Z_A = -L * cos(alpha) - D * sin(alpha);
M_A = c_M * dyn_pres * S * c;

% Thrust forces
X_T = T * cos(a_T);
Z_T = -T * sin(a_T);
M_T = T * z_T;

% State deltas
syms Delta_u Delta_w Delta_q Delta_theta Delta_h real;
% Input deltas
syms Delta_T Delta_delta real;
% Derivatives deltas
syms Delta_w_dot;
% Equilibrium values
syms u_eq w_eq q_eq theta_eq h_eq real;
syms T_eq delta_eq real;

lin = @(symbol, var) subs(diff(symbol, var), [u, w, q, theta, h, T, delta], [u_eq, w_eq, q_eq, theta_eq, h_eq, T_eq, delta_eq]);

X_u = 1 / m * lin(X_A, u);
X_w = 1 / m * lin(X_A, w);
X_delta = 1 / m * lin(X_A, delta);
X_T = 1 / m * lin(X_T, T);
Z_u = 1 / m * lin(Z_A, u);
Z_w = 1 / m * lin(Z_A, w);
Z_q = 1 / m * lin(Z_A, q);
Z_w_dot = 1 / m * lin(Z_A, w_dot);
Z_delta = 1 / m * lin(Z_A, delta);
Z_T = 1 / m * lin(Z_T, T);
M_u = 1 / I_yy * lin(M_A, u);
M_w = 1 / I_yy * lin(M_A, w);
M_q = 1 / I_yy * lin(M_A, q);
M_w_dot = 1 / I_yy * lin(M_A, w_dot);
M_delta = 1 / I_yy * lin(M_A, delta);
M_T = 1 / I_yy * lin(M_T, T);

X_u = simplify(X_u);
X_w = simplify(X_w);
X_delta = simplify(X_delta);
X_T = simplify(X_T);
Z_u = simplify(Z_u);
Z_w = simplify(Z_w);
Z_q = simplify(Z_q);
Z_w_dot = simplify(Z_w_dot);
Z_delta = simplify(Z_delta);
Z_T = simplify(Z_T);
M_u = simplify(M_u);
M_w = simplify(M_w);
M_q = simplify(M_q);
M_w_dot = simplify(M_w_dot);
M_delta = simplify(M_delta);
M_T = simplify(M_T);

disp('X_u = '), disp(X_u);
disp('X_w = '), disp(X_w);
disp('X_delta = '), disp(X_delta);
disp('X_T = '), disp(X_T);
disp('Z_u = '), disp(Z_u);
disp('Z_w = '), disp(Z_w);
disp('Z_q = '), disp(Z_q);
disp('Z_w_dot = '), disp(Z_w_dot);
disp('Z_delta = '), disp(Z_delta);
disp('Z_T = '), disp(Z_T);
disp('M_u = '), disp(M_u);
disp('M_w = '), disp(M_w);
disp('M_q = '), disp(M_q);
disp('M_w_dot = '), disp(M_w_dot);
disp('M_delta = '), disp(M_delta);
disp('M_T = '), disp(M_T);

%% From this point on, replace the linearized expressions with symbols
syms X_u X_w X_delta X_T Z_u Z_w Z_q Z_w_dot Z_delta Z_T M_u M_w M_q M_w_dot M_delta M_T real;

Delta_X_A = m * (X_u * Delta_u + X_w * Delta_w + X_delta * Delta_delta);
Delta_Z_A = m * (Z_u * Delta_u + Z_w * Delta_w + Z_q * Delta_q + Z_w_dot * Delta_w_dot + Z_delta * Delta_delta);
Delta_M_A = I_yy * (M_u * Delta_u + M_w * Delta_w + M_q * Delta_q + M_w_dot * Delta_w_dot + M_delta * Delta_delta);

Delta_X_T = m * (X_T * Delta_T);
Delta_Z_T = m * (Z_T * Delta_T);
Delta_M_T = I_yy * (M_T * Delta_T);

% Linearize the EOM
Delta_u_dot = -q_eq * Delta_w - w_eq * Delta_q - g * cos(theta_eq) * Delta_theta + 1 / m * Delta_X_A + 1 / m * Delta_X_T;
Delta_w_dot = q_eq * Delta_u + u_eq * Delta_q - g * sin(theta_eq) * Delta_theta + 1 / m * Delta_Z_A + 1 / m * Delta_Z_T;
Delta_q_dot = 1 / I_yy * Delta_M_A + 1 / I_yy * Delta_M_T;
Delta_theta_dot = Delta_q;

Delta_u_dot = simplify(Delta_u_dot);
Delta_w_dot = simplify(Delta_w_dot);
Delta_q_dot = simplify(Delta_q_dot);
Delta_theta_dot = simplify(Delta_theta_dot);
