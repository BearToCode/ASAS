function PD_study(params, max_overshoot, max_peak_time)
    % PD_STUDY Study the performance of a PI controller for the inverted pendulum system.
    %
    % PD_study(params, max_overshoot, max_peak_time)
    %
    % Input arguments:
    % params            [Struct]           parameters for the inverted pendulum system
    % - M               [1x1]              mass of the cart                                 [kg]
    % - m               [1x1]              mass of the pendulum                             [kg]
    % - l               [1x1]              length of the pendulum                           [m]
    % - g               [1x1]              gravity acceleration                             [m/s^2]
    % max_overshoot     [1x1]              maximum allowed overshoot                        [%]
    % max_peak_time     [1x1]              maximum allowed peak time                        [s]

    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;

    max_Kp = 1000;
    max_Kd = 500;

    [Kp_space, Kd_space] = meshgrid(0.1:max_Kp, 0.1:max_Kd);
    [Kp_space_lr, Kd_space_lr] = meshgrid(linspace(0.1, max_Kp, 20), linspace(0.1, max_Kd, 20)); % Low resolution for overshoot and peak time

    pole = @(Kp, Kd) (-Kd + sqrt(Kd ^ 2 - 4 * M * l * (Kp - (M + m) * g))) / (2 * M * l);

    re = @(Kp, Kd) real(pole(Kp, Kd));
    im = @(Kp, Kd) imag(pole(Kp, Kd));

    overshoot = @(Kp, Kd) 100 * exp(re(Kp, Kd) * pi / im(Kp, Kd));
    peak_time = @(Kp, Kd) pi / im(Kp, Kd);

    Z_stability = arrayfun(@(Kp, Kd) re(Kp, Kd) < 0, Kp_space, Kd_space);
    Z_overshoot = arrayfun(@(Kp, Kd) overshoot(Kp, Kd) < max_overshoot, Kp_space, Kd_space);
    Z_peak_time = arrayfun(@(Kp, Kd) peak_time(Kp, Kd) < max_peak_time, Kp_space, Kd_space);

    Z = Z_stability + Z_overshoot + Z_peak_time; % Combine inequalities

    figure;
    subplot(2, 2, 1)
    surf(Kp_space_lr, Kd_space_lr, arrayfun(@(Kp, Kd) overshoot(Kp, Kd), Kp_space_lr, Kd_space_lr));
    xlabel('$k_P$', 'Interpreter', 'latex')
    ylabel('$k_D$', 'Interpreter', 'latex')
    zlabel('Overshoot [%]', 'Interpreter', 'latex')
    title('Overshoot', Interpreter = 'latex')
    view(120, 45);

    subplot(2, 2, 2)
    surf(Kp_space_lr, Kd_space_lr, arrayfun(@(Kp, Kd) peak_time(Kp, Kd), Kp_space_lr, Kd_space_lr));
    xlabel('$k_P$', 'Interpreter', 'latex')
    ylabel('$k_D$', 'Interpreter', 'latex')
    zlabel('Peak time [s]', 'Interpreter', 'latex')
    title('Peak time', Interpreter = 'latex')
    view(60, 45);

    subplot(2, 2, [3, 4])
    contourf(Kp_space, Kd_space, Z, 'DisplayName', 'Stability and performance region', 'LineColor', 'none');
    xlabel('$k_P$', 'Interpreter', 'latex')
    ylabel('$k_D$', 'Interpreter', 'latex')
    title('Stability and performance region', Interpreter = 'latex')
end
