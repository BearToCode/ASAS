function [F_tt, F_xx, F_tx, F_xt] = ipend_full_control(G_x, G_t, R_x, R_t)
    % IPEND_FULL_CONTROL computes the full control for the inverted pendulum
    %
    % [F_tt, F_xx, F_tx, F_xt] = ipend_full_control(G_x, G_t, R_x, R_t)
    %
    % Input arguments:
    % G_x         [TransferFunction] transfer function of the cart position
    % G_t         [TransferFunction] transfer function of the pendulum angle
    % R_x         [TransferFunction] transfer function of the cart controller
    % R_t         [TransferFunction] transfer function of the pendulum controller
    %
    % Output arguments:
    % F_tt       [TransferFunction] transfer function theta ref -> theta
    % F_xx       [TransferFunction] transfer function x ref -> x
    % F_tx       [TransferFunction] transfer function theta ref -> x
    % F_xt       [TransferFunction] transfer function x ref -> theta

    F_tt = G_t * R_t / (1 + G_x * R_x + G_t * R_t);
    F_xx = G_x * R_x / (1 + G_x * R_x + G_t * R_t);
    F_tx = G_x * R_t / (1 + G_x * R_x + G_t * R_t);
    F_xt = G_t * R_x / (1 + G_x * R_x + G_t * R_t);

end
