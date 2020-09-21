function [Q, R, P, r] = MyQRP(A)
zero_tol = 1e-15;
[m, n] = size(A);
Q = eye(m);
Pp = (1:n);
P = eye(n);
r = min(m,n);

% Compute the norms (squared) of each column
col_norms = vecnorm(A).^2;

for j = 1:n
    % Find the column with max norm
    [~, maxindex] = max(col_norms(j:end));
    maxindex = maxindex + j - 1;
    
    if (max(abs(A(j:end, maxindex))) <= zero_tol) % Rank trigger
        r = j - 1;
        
        if r < n
           A(r+1:end, r+1:end) = 0; 
        end
        break
        
    end
    
    % Update the permutation, the matrix and the norm vector
    tmpP = Pp(j);
    tmpA = A(:, j);
    tmpN = col_norms(j);
    
    Pp(j) = Pp(maxindex);
    A(:, j) = A(:, maxindex);
    col_norms(j) = col_norms(maxindex);
    
    Pp(maxindex) = tmpP;
    A(:, maxindex) = tmpA;
    col_norms(maxindex) = tmpN;
    
    
    % Compute the householder vector %
    [v, s] = householder_vector(A(j:end, j)); 
    
    % Plug in the first column manually %
    A(j,j) = s; 
    A(j+1:end,j) = 0;
    
    A(j:end,j+1:end) = A(j:end,j+1:end) - 2*v*(v'*A(j:end,j+1:end));
    
    Q(:, j:end) = Q(:, j:end) - Q(:,j:end)*v*2*v';
    
    % Update the column norms by subtracting the squared elements
    col_norms(j+1:end) = col_norms(j+1:end) - (A(j, j+1:end)).^2;
end

if (m == n)
    Q(:, m) = -Q(:, m);
    A(n,n) = -A(n,n);
end

P = P(:,Pp);
R = A;

end