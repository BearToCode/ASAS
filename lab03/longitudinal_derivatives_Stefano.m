function stability = longitudinal_derivatives_Stefano(params, aer, x_eq, input_eq)
    % LONGITUDINAL_DERIVATIVES computes the longitudinal stability derivatives
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

    delta = input_eq(1); 
    %T_eq = input_eq(2);
    alpha = atan(w_eq / u_eq);

    c_L = c_L0 + c_L_alpha * alpha + c_L_delta * delta;
    c_D = c_D0 + c_D_alpha * alpha + c_D_alpha2 * alpha ^ 2;
    c_M = c_M0 + c_M_alpha * alpha + c_M_q * (c / sqrt(u_eq ^ 2 * w_eq ^ 2)) * q_eq + c_M_delta * delta;
   
    stability.X_u = (rho*u_eq*S*(c_L*sin(alpha)+c_D*cos(alpha))-1/2*rho*w_eq*S*(c_L_alpha*sin(alpha)+c_L*cos(alpha)+(c_D_alpha+2*c_D_alpha2*alpha)*cos(alpha)-c_D*sin(alpha)))/m;
    stability.X_w = (rho*w_eq*S*(c_L*sin(alpha)+c_D*cos(alpha))+1/2*rho*u_eq*S*(c_L_alpha*sin(alpha)+c_L*cos(alpha)+(c_D_alpha+2*c_D_alpha2*alpha)*cos(alpha)-c_D*sin(alpha)))/m;
    stability.X_delta = 1/2*rho*(u_eq^2+w_eq^2)*S*c_L_delta*sin(alpha)/m;
    stability.X_T = cos(a_T) / m;
    stability.Z_u = (-rho*u_eq*S*(c_L*cos(alpha)+c_D*sin(alpha))+1/2*rho*w_eq*S*(c_L_alpha*cos(alpha)-c_L*sin(alpha)+(c_D_alpha+2*c_D_alpha2*alpha)*sin(alpha)+c_D*cos(alpha)))/m;
    stability.Z_w = (-rho*w_eq*S*(c_L*cos(alpha)+c_D*sin(alpha))+1/2*rho*u_eq*S*(c_L_alpha*cos(alpha)-c_L*sin(alpha)+(c_D_alpha+2*c_D_alpha2*alpha)*sin(alpha)+c_D*cos(alpha)))/m;
    stability.Z_q = 0;
    stability.Z_w_dot = 0;
    stability.Z_delta = -1/2*rho*(u_eq^2+w_eq^2)*S*c_L_delta*cos(alpha)/m;
    stability.Z_T = -sin(a_T) / m;
    stability.M_u = (rho*u_eq*S*c*c_M-1/2*rho*w_eq*S*c*c_M_alpha)/I_yy;
    stability.M_w = (rho*w_eq*S*c*c_M+1/2*rho*u_eq*S*c*c_M_alpha)/I_yy;
    stability.M_q = 1/2*rho*(u_eq^2+w_eq^2)*S*c^2*c_M_q/I_yy;
    stability.M_w_dot = 0;
    stability.M_delta = 1/2*rho*(u_eq^2+w_eq^2)*S*c*c_M_delta/I_yy;
    stability.M_T = z_T / I_yy;
end
