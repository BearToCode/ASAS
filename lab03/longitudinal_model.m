function f = longitudinal_model(params, aer)
    %% LONGITUDINAL_MODEL Returns the longitudinal model of the aircraft
    %
    % f = longitudinal_model(params, aer)
    %
    % Input arguments:
    % params       [Struct]       parameters for the aircraft model
    % aer          [Struct]       aerodynamic model of the aircraft
    %
    % Output arguments:
    % f           [FunctionHandle] function handle for the longitudinal model

    function x_dot = f_extended(x, input)
        u = x(1);
        w = x(2);
        q = x(3);
        theta = x(4);
        h = x(5);

        delta = input(1);
        T = input(2); % Thrust

        [~, ~, ~, rho, ~, ~] = atmosisa(h);

        if (u ~= 0)
            alpha = atan(w / u);
        else
            alpha = 0;
        end

        V = sqrt(u ^ 2 + w ^ 2);

        c_L = aer.c_L(alpha, delta);
        c_D = aer.c_D(alpha);
        c_M = aer.c_M(alpha, delta, V, q);

        L = 0.5 * rho * V ^ 2 * params.S * c_L;
        D = 0.5 * rho * V ^ 2 * params.S * c_D;

        X_A = L * sin(alpha) - D * cos(alpha);
        Z_A = -L * cos(alpha) - D * sin(alpha);
        M_A = 0.5 * rho * V ^ 2 * params.S * params.c * c_M;

        X_T = T * cos(params.a_T);
        Z_T = -T * sin(params.a_T);
        M_T = T * params.z_T;

        u_dot = -q * w - params.g * sin(theta) + 1 / params.m * X_A + 1 / params.m * X_T;
        w_dot = q * u + params.g * cos(theta) + 1 / params.m * Z_A + 1 / params.m * Z_T;
        q_dot = 1 / params.I_yy * M_A + 1 / params.I_yy * M_T;
        theta_dot = q;
        h_dot = u * sin(theta) - w * cos(theta);

        x_dot = [u_dot; w_dot; q_dot; theta_dot; h_dot];
    end

    f = @f_extended;
end
