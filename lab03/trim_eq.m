function diff = trim_eq(x, V_trim, h_trim, f)
    % TRIM_EQ Trim equations for the aircraft model
    %
    % diff = trim_eq(x, V_trim, h_trim, f)
    %
    % Input arguments:
    % x           [1x3]             state vector                                [T, alpha, delta]
    % V_trim      [1x1]             trim velocity                               [m/s]
    % h_trim      [1x1]             trim altitude                               [m]
    % f           [FunctionHandle]  function handle for the longitudinal model
    %
    % Output arguments:
    % diff        [1x5]             difference between the model state and the trim conditions

    T = x(1); % thrust [N]
    alpha = x(2); % angle of attack [rad]
    delta = x(3); % elevator deflection [rad]

    u = V_trim * cos(alpha); % velocity in x [m/s]
    w = V_trim * sin(alpha); % velocity in y [m/s]
    q = 0; % pitch rate [rad/s]
    theta = alpha; % pitch angle [rad] equivalent to angle of attack
    h = h_trim; % altitude [m]

    diff = f([u; w; q; theta; h], [delta; T]) - [0; 0; 0; 0; 0]; % steady-state condition
end
