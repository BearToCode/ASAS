function [A, B_u, B_d, C, D_u, D_d] = ipend_control_linear(params)
    % IPEND_LINEAR Return the linearized state space model of the inverted pendulum control system.
    %
    % sys = ipend_control_linear(params)
    %
    % Input arguments:
    % params     [Struct]       parameters for the inverted pendulum system
    % - M        [1x1]          mass of the cart                                [kg]
    % - m        [1x1]          mass of the pendulum                            [kg]
    % - l        [1x1]          length of the pendulum                          [m]
    % - g        [1x1]          gravity acceleration                            [m/s^2]
    %
    % Output arguments:
    % A          [4x4]          state matrix of the linearized system
    % B_u        [4x1]          input matrix of the linearized system
    % B_d        [4x1]          disturbance matrix of the linearized system
    % C          [2x4]          output matrix of the linearized system
    % D_u        [2x1]          direct transmission matrix of the linearized system
    % D_d        [2x1]          direct transmission matrix of the linearized system

    M = params.M;
    m = params.m;
    l = params.l;
    g = params.g;

    A = [
         0 1 0 0;
         0 0 (+m * g / M) 0;
         0 0 0 1;
         0 0 ((1 + m / M) * (g / l)) 0
         ];
    B_u = [0; 1 / M; 0; 1 / (M * l)];
    B_d = [0; 1 / M; 0; 1 / (M * l)];
    C = [1 0 0 0;
         0 0 1 0];
    D_u = [0; 0];
    D_d = [0; 0];
end
