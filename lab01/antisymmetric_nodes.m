function [xi] = antisymmetric_nodes(n_modes)
    % ANTISYMMETRIC_NODES returns the first n_modes antisymmetric nodes for the
    % sloshing model.
    %
    % xi = antisymmetric_nodes(n_modes) returns the first n_modes antisymmetric
    %
    % Input arguments:
    % n_modes       [1x1]       number of modes to calculate
    %
    % Output arguments:
    % xi            [nx1]       antisymmetric nodes

    xi = zeros(n_modes, 1);
    xi(1:3) = [1.841; 5.329; 8.531];

    xi = xi(1:n_modes);

    for k = 4:n_modes
        xi(k) = xi(k - 1) + pi;
    end

end
