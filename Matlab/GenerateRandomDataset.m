    %   M       N       Rk      Cnt     k1      k2      k3  ...
Batches = { ...
    genBatch(10, 20, 10, 3, 9, 2, 5); ... % Rectangular matrices
    genBatch(20, 40, 20, 3, 19, 2, 9); ... 
    genBatch(40, 80, 40, 3, 39, 2, 9); ...
    genBatch(80, 160, 80, 3, 79, 2, 9); ...
    genBatch(160, 320, 160, 3, 159, 2, 9); ...
    genBatch(320, 640, 320, 3, 319, 2, 9); ...
    genBatch(10, 10, 10, 3, 9, 2, 5); ... %Square Matrices
    genBatch(20, 20, 20, 3, 19, 2, 9); ...
    genBatch(40, 40, 40, 3, 39, 2, 9); ...
    genBatch(80, 80, 80, 3, 79, 2, 9); ...
    genBatch(160, 160, 160, 3, 159, 2, 9); ...
    genBatch(320, 320, 320, 3, 319, 2, 9); ...
    genBatch(10, 20, 7, 3, 6, 2, 5); ... % Rectangular 3/4-rank matrices
    genBatch(20, 40, 15, 3, 14, 2, 9); ... 
    genBatch(40, 80, 30, 3, 29, 2, 9); ...
    genBatch(80, 160, 60, 3, 59, 2, 9); ...
    genBatch(160, 320, 120, 3, 119, 2, 9); ...
    genBatch(320, 640, 240, 3, 239, 2, 9); ...
    genBatch(10, 10, 7, 3, 6, 2, 5); ... % Square 3/4-rank matrices
    genBatch(20, 20, 15, 3, 14, 2, 9); ... 
    genBatch(40, 40, 30, 3, 29, 2, 9); ...
    genBatch(80, 80, 60, 3, 59, 2, 9); ...
    genBatch(160, 160, 120, 3, 119, 2, 9); ...
    genBatch(320, 320, 240, 3, 239, 2, 9); ...
    };

%{
    e.g. A line containing:
    [   100,    200,    90,     7,      89,     50 ,    10]
    Generates a batch consisting of:
        - 7 Matrices (Cnt)
        - Each having dimension 100 x 200 (M x N)
        - Each having rank 90 (Rk)
        - And asks to give a low rank approximation for each of the ranks:
            - 89 (k1)
            - 50 (k2)
            - 10 (k3)
%}
    
outfolder = 'data_random';

for i = 1:size(Batches, 1)
    CurLine = Batches{i};
    M = CurLine(1);
    N = CurLine(2);
    Rk = CurLine(3);
    Cnt = CurLine(4);
    
    for j = 1:Cnt
        fprintf("Generating batch %d, task %d\n", i, j);
        A = RandRank(M, N, Rk)*1000;
        
        tic % Compute the svd factorization to derive the optimal approx.
        [U,S,V] = svd(A);
        Time_opt = toc;
        
        for h = 5:length(CurLine)
            k = CurLine(h);
            
            Ak_opt = U(:, 1:k)*S(1:k,1:k)*V(:,1:k)';
            Ak_delta = norm(A - Ak_opt, 'fro')/norm(A, 'fro');
            title = sprintf("%s/rnd_%d_%d_%d.mat", outfolder, i, j, (h-4));
            save(title, 'M', 'N', 'Rk', 'k', 'A', 'Ak_opt', 'Ak_delta', 'Time_opt');
        end
    end
end


function [b] = genBatch(M, N, Rk, Cnt, maxrk, minrk, nk)
    step = (maxrk-minrk)/nk;
    step = max(step, 1);
    
    ks = round(minrk:step:maxrk);
    ks(end) = maxrk;
    
    b = [M, N, Rk, Cnt, ks];
end
