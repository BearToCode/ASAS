function [F_den] = ipend_full_control_den(params)
    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;

    F_den = @(Kp_x, Kd_x, Kp_theta, Kd_theta) [- M * l, - Kd_theta - Kd_x * l, - Kp_theta + g * m + M * g - Kp_x * l, + Kd_x * g, + Kp_x * g];
end
