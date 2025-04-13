function [f, g] = ipend_nonlinear(params)
    % IPEND_NONLINEAR Return the nonlinear equations of motion of the inverted pendulum system.
    %
    % [f, g] = ipend_nonlinear(params)
    %
    % Input arguments:
    % params     [Struct]       parameters for the inverted pendulum system
    % - M        [1x1]          mass of the cart                                [kg]
    % - m        [1x1]          mass of the pendulum                            [kg]
    % - l        [1x1]          length of the pendulum                          [m]
    % - g        [1x1]          gravity acceleration                            [m/s^2]
    % - c        [1x1]          damping coefficient of the cart                 [N/m/s]
    % - b        [1x1]          damping coefficient of the pendulum             [Nm/rad/s]
    %
    % Output arguments:
    % f          [Function]     nonlinear EOM of the inverted pendulum system
    %                           f(x, u) = dx/dt
    % g          [Function]     nonlinear EOM of the inverted pendulum system
    %                           g(x, u) = y

    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;
    c = params.c;
    b = params.b;

    f = @(x, u) [
                 x(2); % dx
                 (b / l * x(4) * cos(x(3)) + m * l * x(4) ^ 2 * sin(x(3)) - c * x(2) - m * g * sin(x(3)) * cos(x(3)) + u) / ...
                     (M + m * sin(x(3)) ^ 2); % d2x
                 x(4); % dtheta
                 ((M + m) * g / l * sin(x(3)) - (M / m + 1) * b / l ^ 2 * x(4) - m * l * x(4) ^ 2 * cos(x(3)) * sin(x(3)) + c / l * x(2) * cos(x(3)) - u * cos(x(3)) / l) / ...
                     (M + m * sin(x(3)) ^ 2); % d2theta
                 ];

    g = @(x, u) [
                 x(1); % x
                 x(3); % theta
                 ];
end
