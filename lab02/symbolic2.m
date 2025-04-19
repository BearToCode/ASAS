clc; clear;

syms t real;
syms m M l g u d real;
syms x(t) theta(t);

dx = diff(x, t);
dtheta = diff(theta, t);

T = 1/2 * (m + M) * dx ^ 2 +1/2 * m * l ^ 2 * dtheta ^ 2 - m * l * dx * dtheta * cos(theta);
V = m * g * l * cos(theta);

Q_x = u + d;
Q_theta = 0;

f_x = diff(diff(T, dx), t) - diff(T, x) + diff(V, x) - Q_x;
f_theta = diff(diff(T, dtheta), t) - diff(T, theta) + diff(V, theta) - Q_theta;

eq_x = f_x == 0;
eq_theta = f_theta == 0;

sys = simplify([eq_x; eq_theta]);

syms d2x d2theta real;
sys = subs(sys, [diff(x, t, t), diff(theta, t, t)], [d2x, d2theta]);

solved_sys = solve(sys, [d2x, d2theta]);
disp("Second-order EOM:");
latex(simplify(solved_sys.d2x))
latex(simplify(solved_sys.d2theta))

% Linearization

x0 = 0; % initial cart position
theta0 = 0; % initial pendulum angle
dx0 = 0; % initial cart velocity
dtheta0 = 0; % initial pendulum angular velocity

syms delta_x delta_theta delta_dx delta_dtheta delta_d2x delta_d2theta real;

f_x_eq = eval_at_equilibrium(f_x, x, theta, t, x0, theta0, dx0, dtheta0);
f_x_x = eval_at_equilibrium(diff(f_x, x), x, theta, t, x0, theta0, dx0, dtheta0);
f_x_dx = eval_at_equilibrium(diff(f_x, diff(x, t)), x, theta, t, x0, theta0, dx0, dtheta0);
f_x_d2x = eval_at_equilibrium(diff(f_x, diff(x, t, t)), x, theta, t, x0, theta0, dx0, dtheta0);
f_x_theta = eval_at_equilibrium(diff(f_x, theta), x, theta, t, x0, theta0, dx0, dtheta0);
f_x_dtheta = eval_at_equilibrium(diff(f_x, diff(theta, t)), x, theta, t, x0, theta0, dx0, dtheta0);
f_x_d2theta = eval_at_equilibrium(diff(f_x, diff(theta, t, t)), x, theta, t, x0, theta0, dx0, dtheta0);

f_theta_eq = eval_at_equilibrium(f_theta, x, theta, t, x0, theta0, dx0, dtheta0);
f_theta_x = eval_at_equilibrium(diff(f_theta, x), x, theta, t, x0, theta0, dx0, dtheta0);
f_theta_dx = eval_at_equilibrium(diff(f_theta, diff(x, t)), x, theta, t, x0, theta0, dx0, dtheta0);
f_theta_d2x = eval_at_equilibrium(diff(f_theta, diff(x, t, t)), x, theta, t, x0, theta0, dx0, dtheta0);
f_theta_theta = eval_at_equilibrium(diff(f_theta, theta), x, theta, t, x0, theta0, dx0, dtheta0);
f_theta_dtheta = eval_at_equilibrium(diff(f_theta, diff(theta, t)), x, theta, t, x0, theta0, dx0, dtheta0);
f_theta_d2theta = eval_at_equilibrium(diff(f_theta, diff(theta, t, t)), x, theta, t, x0, theta0, dx0, dtheta0);

f_x_lin = f_x_eq + f_x_x * delta_x + f_x_dx * delta_dx + f_x_d2x * delta_d2x + ...
    f_x_theta * delta_theta + f_x_dtheta * delta_dtheta + f_x_d2theta * delta_d2theta;

f_theta_lin = f_theta_eq + f_theta_x * delta_x + f_theta_dx * delta_dx + f_theta_d2x * delta_d2x + ...
    f_theta_theta * delta_theta + f_theta_dtheta * delta_dtheta + f_theta_d2theta * delta_d2theta;

lin_sys = [f_x_lin; f_theta_lin];

solved_lin_sys = solve(lin_sys, [delta_d2x, delta_d2theta]);
disp("Linearized EOM:");
latex(simplify(solved_lin_sys.delta_d2x))
latex(simplify(solved_lin_sys.delta_d2theta))

function static = eval_at_equilibrium(sym, x_var, theta_var, t_var, x0, theta0, dx0, dtheta0)
    static = subs(sym, {x_var, theta_var, diff(x_var, t_var), diff(theta_var, t_var)}, {x0, theta0, dx0, dtheta0});
    static = subs(static, {diff(x_var, t_var, t_var), diff(theta_var, t_var, t_var)}, {0, 0});
    static = simplify(static);
end

syms s;

A = [
     0 1 0 0;
     0 0 (+m * g / M) 0;
     0 0 0 1;
     0 0 ((1 + m / M) * (g / l)) 0
     ];
B = [0; 1 / M; 0; 1 / (M * l)];
C = [1 0 0 0;
     0 0 1 0];
D = [0; 0];

G_x = C(1, :) * inv(s * eye(4) - A) * B + D(1, :);
G_theta = C(2, :) * inv(s * eye(4) - A) * B + D(2, :);

G_x = simplify(G_x);
G_theta = simplify(G_theta);

syms Kp_x Kd_x Kp_theta Kd_theta real;
R_x = Kp_x + Kd_x * s;
R_theta = Kp_theta + Kd_theta * s;

[F_tt, F_xx, F_tx, F_xt] = ipend_full_control(G_x, G_theta, R_x, R_theta);
[~, F_tt_den] = numden(simplify(F_tt));
[~, F_xx_den] = numden(simplify(F_xx));
[~, F_tx_den] = numden(simplify(F_tx));
[~, F_xt_den] = numden(simplify(F_xt));

% AlL THE SAME!!
