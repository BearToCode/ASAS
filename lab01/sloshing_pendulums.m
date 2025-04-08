function pendulums = sloshing_pendulums(params, n_modes)
    % SLOSHING_PENDULUMS calculates the parameters for the sloshing pendulums
    % model.
    %
    % pendulums = sloshing_pendulums(params, n_modes) returns the parameters for
    % the sloshing pendulums model.
    %
    % Input arguments:
    % params        [Struct]    parameters for the model
    % - h           [1x1]       height of the tank                          [m]
    % - d           [1x1]       diameter of the tank                        [m]
    % - g           [1x1]       gravity acceleration                        [m/s^2]
    % - density     [1x1]       density of the fluid                        [kg/m^3]
    % n_modes       [1x1]       number of modes to calculate                [-]
    %
    % Output arguments:
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

    M = params.d ^ 2 * pi / 4 * params.h * params.density;

    xi = antisymmetric_nodes(n_modes);

    % Pendulum masses
    m = M * params.d * tanh(2 * xi * params.h / params.d) ./ xi ./ (xi .^ 2 - 1) / params.h;
    m0 = M - sum(m);

    % Pendulums lengths
    L = params.d ./ (2 * xi .* tanh(2 * xi * params.h / params.d));

    % Pendulum heights, taken from center of gravity
    H = params.h / 2 - params.d / 2 ./ xi .* (tanh(xi * params.h / params.d) - 1 ./ sinh(2 * xi * params.h / params.d));
    H0 = sum(m .* (H - L)) / m0;

    % Natural frequencies
    w_n = sqrt(params.g ./ L);
    f_n = w_n ./ (2 * pi);

    pendulums.n = n_modes;
    pendulums.xi = xi;
    pendulums.m = m;
    pendulums.m0 = m0;
    pendulums.L = L;
    pendulums.H = H;
    pendulums.H0 = H0;
    pendulums.w_n = w_n;
    pendulums.f_n = f_n;
end
