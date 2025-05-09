function model = aerodynamic_model(params)
    % AERODYNAMIC_MODEL Returns the aerodynamic model of the aircraft
    %
    % model = aerodynamic_model(params)
    %
    % Input arguments:
    % params     [Struct]       parameters for the aircraft model
    % - m        [1x1]          mass of the aircraft                            [Kg]
    % - I_yy     [1x1]          moment of inertia around the y-axis             [Kg*m^2]
    % - S        [1x1]          wing area                                       [m^2]
    % - c        [1x1]          mean aerodynamic chord                          [m]
    % - a_T      [1x1]          thrustline angle                                [rad]
    % - z_T      [1x1]          thrustline vertical distance                    [m]
    %
    % Output arguments:
    % model         [Struct]        aerodynamic model of the aircraft
    % - c_L0        [1x1]           lift coefficient at zero angle of attack         [-]
    % - c_L_alpha   [1x1]           lift coefficient due to angle of attack          [rad^-1]
    % - c_L_delta   [1x1]           lift coefficient due to elevator deflection      [rad^-1]
    % - c_D0        [1x1]           drag coefficient at zero angle of attack         [-]
    % - c_D_alpha   [1x1]           drag coefficient due to angle of attack          [rad^-1]
    % - c_D_alpha2  [1x1]           drag coefficient due to angle of attack^2        [rad^-2]
    % - c_M0        [1x1]           moment coefficient at zero angle of attack       [-]
    % - c_M_alpha   [1x1]           moment coefficient due to angle of attack        [rad^-1]
    % - c_Mq        [1x1]           moment coefficient due to pitch rate             [rad^-1]
    % - c_M_delta   [1x1]           moment coefficient due to elevator deflection    [rad^-1]
    % - c_L         [function]      lift coefficient function                        [rad^-1]
    % - c_D         [function]      drag coefficient function                        [rad^-1]
    % - c_M         [function]      moment coefficient function                      [rad^-1]

    model.c_L0 = 0.895;
    model.c_L_alpha = 5.01;
    model.c_L_delta = 0.722;
    model.c_D0 = 0.177;
    model.c_D_alpha = 0.232;
    model.c_D_alpha2 = 1.393;
    model.c_M0 = -0.046;
    model.c_M_alpha = -1.087;
    model.c_M_q = -7.055;
    model.c_M_delta = -1.88;

    model.c_L = @(alpha, delta) model.c_L0 + model.c_L_alpha * alpha + model.c_L_delta * delta;
    model.c_D = @(alpha) model.c_D0 + model.c_D_alpha * alpha + model.c_D_alpha2 * alpha ^ 2;
    model.c_M = @(alpha, delta, V, q) model.c_M0 + model.c_M_alpha * alpha + model.c_M_q * (params.c / V) * q + model.c_M_delta * delta;
end
