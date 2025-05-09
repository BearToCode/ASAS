function stability = longitudinal_derivatives(params, aer, x_eq, input_eq)
    %% LONGITUDINAL_DERIVATIVES computes the longitudinal stability derivatives
    %
    % stability = longitudinal_derivatives(params, aer, x_eq, input_eq)
    %
    % Input arguments:
    % params     [Struct]       parameters for the aircraft model
    % aer        [Struct]       aerodynamic model of the aircraft
    % x_eq       [1x5]          equilibrium state vector
    %
    % Output arguments:
    % stability         [Struct]       longitudinal stability derivatives
    % - X_u             [1x1]          longitudinal force derivative w.r.t. u         [N/(m/s)]
    % - X_w             [1x1]          longitudinal force derivative w.r.t. w         [N/(m/s)]
    % - X_delta         [1x1]          longitudinal force derivative w.r.t. delta     [N/(rad)]
    % - X_T             [1x1]          longitudinal force derivative w.r.t. T         [N/(N)]
    % - Z_u             [1x1]          longitudinal force derivative w.r.t. u         [N/(m/s)]
    % - Z_w             [1x1]          longitudinal force derivative w.r.t. w         [N/(m/s)]
    % - Z_q             [1x1]          longitudinal force derivative w.r.t. q         [N/(rad/s)]
    % - Z_delta         [1x1]          longitudinal force derivative w.r.t. delta     [N/(rad)]
    % - Z_T             [1x1]          longitudinal force derivative w.r.t. T         [N/(N)]
    % - M_u             [1x1]          longitudinal moment derivative w.r.t. u        [Nm/(m/s)]
    % - M_w             [1x1]          longitudinal moment derivative w.r.t. w        [Nm/(m/s)]
    % - M_q             [1x1]          longitudinal moment derivative w.r.t. q        [Nm/(rad/s)]
    % - M_delta         [1x1]          longitudinal moment derivative w.r.t. delta    [Nm/(rad)]
    % - M_T             [1x1]          longitudinal moment derivative w.r.t. T        [Nm/(N)]
    % - M_w_dot         [1x1]          longitudinal moment derivative w.r.t. w_dot    [Nm/(m/s^2)]

    m = params.m;
    I_yy = params.I_yy;
    S = params.S;
    c = params.c;
    a_T = params.a_T;
    z_T = params.z_T;
    % g = params.g;

    c_L0 = aer.c_L0;
    c_L_alpha = aer.c_L_alpha;
    c_L_delta = aer.c_L_delta;
    c_D0 = aer.c_D0;
    c_D_alpha = aer.c_D_alpha;
    c_D_alpha2 = aer.c_D_alpha2;
    c_M0 = aer.c_M0;
    c_M_alpha = aer.c_M_alpha;
    c_M_q = aer.c_M_q;
    c_M_delta = aer.c_M_delta;

    u_eq = x_eq(1);
    w_eq = x_eq(2);
    q_eq = x_eq(3);
    % theta_eq = x_eq(4);
    h_eq = x_eq(5);

    [~, ~, ~, rho, ~, ~] = atmosisa(h_eq);

    delta_eq = input_eq(1);
    % T_eq = input_eq(2);

    stability.X_u =- (S * rho * (u_eq ^ 2) ^ (3/2) * (2 * c_D0 * u_eq ^ 2 + c_D0 * w_eq ^ 2 + c_L_alpha * w_eq ^ 2 + 2 * c_D_alpha2 * u_eq ^ 2 * atan(w_eq / u_eq) ^ 2 + c_D_alpha2 * w_eq ^ 2 * atan(w_eq / u_eq) ^ 2 + 2 * c_D_alpha * u_eq ^ 2 * atan(w_eq / u_eq) + c_D_alpha * w_eq ^ 2 * atan(w_eq / u_eq) - c_D_alpha * u_eq * w_eq - c_L0 * u_eq * w_eq - c_L_delta * delta_eq * u_eq * w_eq - 2 * c_D_alpha2 * u_eq * w_eq * atan(w_eq / u_eq) - c_L_alpha * u_eq * w_eq * atan(w_eq / u_eq))) / (2 * m * u_eq ^ 3 * (u_eq ^ 2 + w_eq ^ 2) ^ (1/2));
    stability.X_w = (S * rho * (u_eq ^ 2) ^ (3/2) * (c_L0 * u_eq ^ 2 - c_D_alpha * u_eq ^ 2 + 2 * c_L0 * w_eq ^ 2 - 2 * c_D_alpha2 * u_eq ^ 2 * atan(w_eq / u_eq) + c_L_alpha * u_eq ^ 2 * atan(w_eq / u_eq) + 2 * c_L_alpha * w_eq ^ 2 * atan(w_eq / u_eq) - c_D0 * u_eq * w_eq + c_L_alpha * u_eq * w_eq + c_L_delta * delta_eq * u_eq ^ 2 + 2 * c_L_delta * delta_eq * w_eq ^ 2 - c_D_alpha2 * u_eq * w_eq * atan(w_eq / u_eq) ^ 2 - c_D_alpha * u_eq * w_eq * atan(w_eq / u_eq))) / (2 * m * u_eq ^ 3 * (u_eq ^ 2 + w_eq ^ 2) ^ (1/2));
    stability.X_delta = (S * c_L_delta * rho * w_eq * abs(u_eq) * (u_eq ^ 2 + w_eq ^ 2) ^ (1/2)) / (2 * m * u_eq);
    stability.X_T = cos(a_T) / m;
    stability.Z_u =- (S * rho * (u_eq ^ 2) ^ (3/2) * (2 * c_L0 * u_eq ^ 2 - c_D_alpha * w_eq ^ 2 + c_L0 * w_eq ^ 2 + 2 * c_L_alpha * u_eq ^ 2 * atan(w_eq / u_eq) - 2 * c_D_alpha2 * w_eq ^ 2 * atan(w_eq / u_eq) + c_L_alpha * w_eq ^ 2 * atan(w_eq / u_eq) + c_D0 * u_eq * w_eq - c_L_alpha * u_eq * w_eq + 2 * c_L_delta * delta_eq * u_eq ^ 2 + c_L_delta * delta_eq * w_eq ^ 2 + c_D_alpha2 * u_eq * w_eq * atan(w_eq / u_eq) ^ 2 + c_D_alpha * u_eq * w_eq * atan(w_eq / u_eq))) / (2 * m * u_eq ^ 3 * (u_eq ^ 2 + w_eq ^ 2) ^ (1/2));
    stability.Z_w =- (S * rho * (u_eq ^ 2) ^ (3/2) * (c_D0 * u_eq ^ 2 + c_L_alpha * u_eq ^ 2 + 2 * c_D0 * w_eq ^ 2 + c_D_alpha2 * u_eq ^ 2 * atan(w_eq / u_eq) ^ 2 + 2 * c_D_alpha2 * w_eq ^ 2 * atan(w_eq / u_eq) ^ 2 + c_D_alpha * u_eq ^ 2 * atan(w_eq / u_eq) + 2 * c_D_alpha * w_eq ^ 2 * atan(w_eq / u_eq) + c_D_alpha * u_eq * w_eq + c_L0 * u_eq * w_eq + c_L_delta * delta_eq * u_eq * w_eq + 2 * c_D_alpha2 * u_eq * w_eq * atan(w_eq / u_eq) + c_L_alpha * u_eq * w_eq * atan(w_eq / u_eq))) / (2 * m * u_eq ^ 3 * (u_eq ^ 2 + w_eq ^ 2) ^ (1/2));
    stability.Z_q = 0;
    stability.Z_w_dot = 0;
    stability.Z_delta =- (S * c_L_delta * rho * abs(u_eq) * (u_eq ^ 2 + w_eq ^ 2) ^ (1/2)) / (2 * m);
    stability.Z_T = -sin(a_T) / m;
    stability.M_u =- ((S * c * rho * (c_M_alpha * w_eq * abs(w_eq) * (u_eq ^ 2) ^ (3/2) + c * c_M_q * q_eq * u_eq ^ 3 + c * c_M_q * q_eq * u_eq * w_eq ^ 2)) / (2 * abs(w_eq) * (u_eq ^ 2) ^ (3/2)) - (S * c * rho * u_eq * (c_M0 * abs(u_eq) * abs(w_eq) + c * c_M_q * q_eq + c_M_alpha * atan(w_eq / u_eq) * abs(u_eq) * abs(w_eq) + c_M_delta * delta_eq * abs(u_eq) * abs(w_eq))) / (abs(u_eq) * abs(w_eq))) / I_yy;
    stability.M_w = ((S * c * rho * (u_eq ^ 2 + w_eq ^ 2) * (c_M_alpha / (u_eq * (w_eq ^ 2 / u_eq ^ 2 + 1)) - (c * c_M_q * q_eq * u_eq ^ 2 * w_eq) / (u_eq ^ 2 * w_eq ^ 2) ^ (3/2))) / 2 + S * c * rho * w_eq * (c_M0 + c_M_delta * delta_eq + c_M_alpha * atan(w_eq / u_eq) + (c * c_M_q * q_eq) / (u_eq ^ 2 * w_eq ^ 2) ^ (1/2))) / I_yy;
    stability.M_q = (S * c ^ 2 * c_M_q * rho * (u_eq ^ 2 + w_eq ^ 2)) / (2 * I_yy * abs(u_eq) * abs(w_eq));
    stability.M_w_dot = 0;
    stability.M_delta = (S * c * c_M_delta * rho * (u_eq ^ 2 + w_eq ^ 2)) / (2 * I_yy);
    stability.M_T = z_T / I_yy;
end
