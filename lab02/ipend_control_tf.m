function [G_x, G_theta] = ipend_control_tf(params)
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
    %
    % Output arguments:
    % G_x        [TransferFunction] transfer function of the cart position
    % G_theta    [TransferFunction] transfer function of the pendulum angle

    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;

    s = tf('s');

    G_x = (- l * s ^ 2 + g) / (s ^ 2 * (- M * l * s ^ 2 + M * g + g * m));
    G_theta = -1 / (- M * l * s ^ 2 + M * g + g * m);
end
