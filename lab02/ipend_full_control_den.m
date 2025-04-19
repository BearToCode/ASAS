function [F_den] = ipend_full_control_den(params)
    % IPEND_FULL_CONTROL_DEN return the denominator of the transfer function for the full control of the inverted pendulum on a cart.
    %
    % This function computes the denominator of the transfer function for the full control of the inverted pendulum on a cart.
    %
    % Input arguments:
    % params        [struct]        Structure containing the parameters of the inverted pendulum on a cart:
    %
    % Output arguments:
    % F_den         [Function]      Returns the polynomial of the denominator of the transfer function

    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;

    F_den = @(Kp_x, Kd_x, Kp_theta, Kd_theta) [- M * l, - Kd_theta - Kd_x * l, - Kp_theta + g * m + M * g - Kp_x * l, + Kd_x * g, + Kp_x * g];
end
