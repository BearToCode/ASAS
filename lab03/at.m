function x = at(A, x, y)
    % Basic function that allows indexing into a matrix. The reason for itx
    % existence is MATLAB not allowing concatenation of indices, for example
    % when you want to index the returned value of a function.
    %
    % x = at(A, x, y)
    %
    % Input arguments:
    % - A: matrix
    % - x: row vector of indices
    % - y: column vector of indices

    x = A(x, y);
end
