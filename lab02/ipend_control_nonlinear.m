function [f, g] = ipend_control_nonlinear(params)
    % IPEND_NONLINEAR Return the nonlinear equations of motion of the inverted control pendulum system.
    %
    % [f, g] = ipend_control_nonlinear(params)
    %
    % Input arguments:
    % params     [Struct]       parameters for the inverted pendulum system
    % - M        [1x1]          mass of the cart                                [kg]
    % - m        [1x1]          mass of the pendulum                            [kg]
    % - l        [1x1]          length of the pendulum                          [m]
    % - g        [1x1]          gravity acceleration                            [m/s^2]
    %
    % Output arguments:
    % f          [Function]     nonlinear EOM of the inverted pendulum control system
    %                           f(x, u, d) = dx/dt
    % g          [Function]     nonlinear EOM of the inverted pendulum control system
    %                           g(x, u, d) = y

    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;

    f = @(x, u, d) [
                    x(2); % dx
                    (m * l * x(4) ^ 2 * sin(x(3)) + m * g * sin(x(3)) * cos(x(3)) + d + u) / ...
                        (M + m * sin(x(3)) ^ 2); % d2x
                    x(4); % dtheta
                    ((M + m) * g / l * sin(x(3)) - m * x(4) ^ 2 * cos(x(3)) * sin(x(3)) + (d + u) * cos(x(3)) / l) / ...
                        (M + m * sin(x(3)) ^ 2); % d2theta
                    ];
    g = @(x, u) [
                 x(1); % x
                 x(3); % theta
                 ];
end
