function [ng] = GradientNorm(A, X, U, V)
[m,n] = size(A);
[~,k] = size(U);

D = A - X;
g = zeros(m*k + k*n, 1);

for a = 1:m
    si = (a-1)*k + 1;
    ei = a*k;
    
    g(si:ei) = V*(D(a,:))';
end

ss = m*k;

for b = 1:n
    si = ss + (b-1)* k + 1;
    ei = ss + b*k;
    
    g(si:ei) = U' * D(:,b);
end

g = -2*g;
ng = norm(g, 'fro');
end

