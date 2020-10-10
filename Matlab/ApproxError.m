% Computes the approximation error using the Frobenius norm for a target
% matrix A, w.r.t. a candidate matrix X

function [n] = ApproxError(A, X)
    d = A - X;
    n = norm(d, 'fro');
end