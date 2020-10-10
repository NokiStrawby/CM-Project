% Utility script to generate random datasets of sparse matrices

%   M       N      Cnt     nk  ...

Prelim_Sparse_Batches = { ...
    [10, 10, 3, 5]; ... %Square Matrices
    [20, 20, 3, 5]; ...
    [40, 40, 3, 5]; ...
    [80, 80, 3, 5]; ...
    [160, 160, 3, 5]; ...
    [320, 320, 3, 5]; ...
    };

Batches = Prelim_Sparse_Batches;
density = 0.60;
    
outfolder = 'data_random';

for i = 1:size(Batches, 1)
    CurLine = Batches{i};
    M = CurLine(1);
    N = CurLine(2);
    Cnt = CurLine(3);
    nk = CurLine(4);
    
    for j = 1:Cnt
        fprintf("Generating batch %d, task %d\n", i, j);
        
        % Generate a matrix of given density
        A = full(sprand(M,N,density))*1000;
        Rk = rank(A);
        
        tic % Compute the svd factorization to derive the optimal approx.
        [U,S,V] = svd(A);
        Time_opt = toc;
        
        ks = equalSpaced(Rk-1, 2, nk);
        
        for h = 1:length(ks)
            k = ks(h);
            Ak_opt = U(:, 1:k)*S(1:k,1:k)*V(:,1:k)';
            Ak_delta = norm(A - Ak_opt, 'fro')/norm(A, 'fro');
            title = sprintf("%s/srnd_%d_%d_%d.mat", outfolder, i, j, h);
            save(title, 'M', 'N', 'Rk', 'k', 'A', 'Ak_opt', 'Ak_delta', 'Time_opt');
        end
    end
end

function [b] = equalSpaced(maxrk, minrk, nk)
    step = (maxrk-minrk)/nk;
    step = max(step, 1);
    
    ks = round(minrk:step:maxrk);
    ks(end) = maxrk;
    b = ks;
end
