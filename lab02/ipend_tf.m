function [G_x, G_theta] = ipend_tf(params)
    % IPEND_TF Return the transfer functions of the inverted pendulum system.
    %
    % [G_x, G_theta] = ipend_tf(params)
    %
    % Input arguments:
    % params     [Struct]           parameters for the inverted pendulum system
    % - M        [1x1]              mass of the cart                                [kg]
    % - m        [1x1]              mass of the pendulum                            [kg]
    % - l        [1x1]              length of the pendulum                          [m]
    % - g        [1x1]              gravity acceleration                            [m/s^2]
    % - c        [1x1]              damping coefficient of the cart                 [N/m/s]
    % - b        [1x1]              damping coefficient of the pendulum             [Nm/rad/s]
    %
    % Output arguments:
    % G_x        [TransferFunction] transfer function of the cart position
    % G_theta    [TransferFunction] transfer function of the pendulum angle

    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;
    c = params.c;
    b = params.b;

    s = tf('s');

    G_x = (m * l ^ 2 * s ^ 2 + g * m * l + b * s) / (s * (M * l ^ 2 * m * s ^ 3 + c * l ^ 2 * m * s ^ 2 + g * l * m ^ 2 * s + M * g * l * m * s + c * g * l * m + b * m * s ^ 2 + M * b * s ^ 2 + b * c * s));
    G_theta = (l * m * s) / (M * l ^ 2 * m * s ^ 3 + c * l ^ 2 * m * s ^ 2 + g * l * m ^ 2 * s + M * g * l * m * s + c * g * l * m + b * m * s ^ 2 + M * b * s ^ 2 + b * c * s);
end
