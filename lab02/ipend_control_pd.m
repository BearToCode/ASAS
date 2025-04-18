function [Kp, Kd, t_a, s_p, p_m] = ipend_control_pd(G)
    [~, G_den] = tfdata(G);
    % Add a minus as matlab returns the numerator as -1
    G_den = -cell2mat(G_den);

    x0 = [1; 1]; % initial guess for Kp and Kd

    gs = GlobalSearch();
    opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
    problem = createOptimProblem('fmincon', ...
        'x0', x0, 'objective', @f, ...
        'lb', [0.01; 0.01], ...
        'options', opts);
    x = run(gs, problem);

    Kp = x(1);
    Kd = x(2);
    [t_a, s_p, p_m] = eval_pd(Kp, Kd, G_den);

    function cost = f(x)
        Kp = x(1);
        Kd = x(2);
        [t_a, s_p, p_m] = eval_pd(Kp, Kd, G_den);
        cost = 0;

        % t_a must be less than 1s
        if t_a > 1
            cost = cost + 1000 * (t_a - 1) ^ 2;
        end

        % s_p must be less than 20%
        if s_p > 20
            cost = cost + 1000 * (s_p - 20) ^ 2;
        end

        % the greater the phase margin the better
        cost = cost - p_m * 1000;

    end

end

function [t_a, s_p, p_m] = eval_pd(Kp, Kd, G_den)
    R_den = [0 Kd Kp];
    L_den = G_den + R_den;

    poles = roots(L_den);
    poles_cc = poles(imag(poles) > 0);

    if isempty(poles_cc)
        t_a = 1e6;
        s_p = 0;
        p_m = 0;
        return;
    end

    pole = poles_cc(1); % Just two poles complex conjugate
    w_n = sqrt(pole * (real(pole) - 1i * imag(pole)));
    e =- (real(pole / w_n));

    t_a = 5 / (w_n * e);
    s_p = 100 * exp(real(pole) / imag(pole) * pi);
    p_m = 100 * e;
end
