function [n] = ApproxError(A, X)
    d = A - X;
    n = norm(d, 'fro');
end