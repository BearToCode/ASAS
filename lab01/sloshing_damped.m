function sys = sloshing_damped(pendulums, params, damping)
    % SLOSHING_DAMPED calculates the state space model for the sloshing
    % pendulums model with damping.
    %
    % sys = sloshing_undamped(pendulums, params, damping) returns the state space model
    % for the sloshing pendulums model with damping.
    %
    % Input arguments:
    % pendulums     [Struct]    parameters for the sloshing pendulums model
    % - n           [1x1]       number of modes                             [n_modes]
    % - xi          [nx1]       antisymmetric nodes                         [-]
    % - m           [nx1]       pendulum masses                             [kg]
    % - m0          [1x1]       mass of the tank                            [kg]
    % - L           [nx1]       pendulum lengths                            [m]
    % - H           [nx1]       pendulum heights                            [m]
    % - H0          [1x1]       height of the tank center of gravity        [m]
    % - w_n         [nx1]       natural frequencies                         [rad/s]
    % - f_n         [nx1]       natural frequencies                         [Hz]
    % params        [Struct]    parameters for the model
    % - h           [1x1]       height of the tank                          [m]
    % - d           [1x1]       diameter of the tank                        [m]
    % - g           [1x1]       gravity acceleration                        [m/s^2]
    % - density     [1x1]       density of the fluid                        [kg/m^3]
    % damping       [1x1]       damping coefficient                         [kg/s]
    %
    % Output arguments:
    % sys           [StateSpace] state space model for the sloshing pendulums
    %                            model without damping

    A = [];

    for i = 1:pendulums.n
        A = blkdiag(A, [0 1; -pendulums.w_n(i) ^ 2, -2 * damping * pendulums.w_n(i)]);
    end

    % Add the final d/dt x0'
    A = blkdiag(A, 0);

    B = zeros(2 * pendulums.n + 1, 1);

    for i = 1:pendulums.n
        B(2 * i - 1:2 * i) = [0; -1 / pendulums.L(i)];
    end

    % Add the final contribution to d/dt x0' from u
    B(2 * pendulums.n + 1) = 1;

    C = zeros(1, 2 * pendulums.n + 1);

    for i = 1:pendulums.n
        C(2 * i - 1) = params.g * pendulums.m(i);
        C(2 * i) = 2 * damping * pendulums.w_n(i) * pendulums.L(i) * pendulums.m(i);
    end

    D = sum(pendulums.m) + pendulums.m0;

    sys = ss(A, B, C, D);
end
