% Generates a random matrix of a given shape and rank, by generating a
% random SVD and truncating it

function X = RandRank(m,n,r)

A = orth(randn(m,m));
B = orth(randn(n,n));
D = diag(rand(r, 1));
X = A(:, 1:r)*D*B(1:r,:);
end