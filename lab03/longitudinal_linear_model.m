function [A, B, C, D] = longitudinal_linear_model(params, stability, x_eq)
    % LONGITUDINAL_LINEAR_MODEL Returns the longitudinal linear model of the aircraft
    %
    % [A, B] = longitudinal_linear_model(params, stability, x_eq)
    %
    % Input arguments:
    % params       [Struct]       parameters for the aircraft model
    % stability    [Struct]       longitudinal stability derivatives
    % x_eq         [1x5]          equilibrium state vector
    %
    % Output arguments:
    % A           [4x4]          state matrix of the longitudinal linear model
    % B           [4x2]          input matrix of the longitudinal linear model
    % C           [6x4]          output matrix of the longitudinal linear model
    % D           [6x2]          feed-forward matrix of the longitudinal linear model

    X_u = stability.X_u;
    X_w = stability.X_w;
    X_delta = stability.X_delta;
    X_T = stability.X_T;
    Z_u = stability.Z_u;
    Z_w = stability.Z_w;
    Z_q = stability.Z_q;
    Z_delta = stability.Z_delta;
    Z_T = stability.Z_T;
    M_u = stability.M_u;
    M_w = stability.M_w;
    M_q = stability.M_q;
    M_delta = stability.M_delta;
    M_T = stability.M_T;

    %     m = params.m;
    %     I_yy = params.I_yy;
    g = params.g;

    u_eq = x_eq(1);
    w_eq = x_eq(2);
    q_eq = x_eq(3);
    theta_eq = x_eq(4);

    A = [
         X_u, (X_w - q_eq), -w_eq, -g * cos(theta_eq);
         (Z_u + q_eq), Z_w, (Z_q + u_eq), -g * sin(theta_eq);
         M_u, M_w, M_q, 0;
         0, 0, 1, 0;
         ];
    B = [
         X_delta, X_T;
         Z_delta, Z_T;
         M_delta, M_T;
         0, 0;
         ];
    C = [
         1, 0, 0, 0;
         0, 1, 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1;
         sin(theta_eq), -cos(theta_eq), 0, (u_eq * cos(theta_eq) + w_eq * sin(theta_eq));
         (-w_eq / u_eq ^ 2) / (1 + (w_eq / u_eq) ^ 2), (1 / u_eq) / (1 + (w_eq / u_eq) ^ 2), 0, 0;
         ];
    D = [0, 0;
         0, 0;
         0, 0;
         0, 0;
         0, 0;
         0, 0];
end
