function sys = ipend_linear(params)
    % IPEND_LINEAR Return the linearized state space model of the inverted pendulum system.
    %
    % sys = ipend_linear(params)
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
    % sys           [StateSpace] state space model for the sloshing pendulums
    %                            model without damping

    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;
    c = params.c;
    b = params.b;

    A = [
         0 1 0 0;
         0 (-c / M) (-m * g / M) (-b / (M * l));
         0 0 0 1;
         0 (- c / (M * l)) (- (1 + m / M) * (g / l)) (- (M + m) / m * b / (M * l ^ 2))
         ];
    B = [0; 1 / M; 0; 1 / (M * l)];
    C = [1 0 0 0;
         0 0 1 0];
    D = [0; 0];

    sys = ss(A, B, C, D);
end
