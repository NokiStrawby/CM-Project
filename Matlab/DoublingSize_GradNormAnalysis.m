% Generates some random square matrices of doubling size,
% executes the algorithm on them and saves the execution
% stats. The execution of this code provides the data 
% for the analysis related to the correlation between 
% the matrix size and the gradient norm flattening.

size_it = 10;
n_task = 1;
max_size = 320;
logsFolder = "Grad_Norm_Analysis";

if ~exist(logsFolder, 'dir')
       mkdir(logsFolder)
end

while size_it <= max_size
    A = randn(size_it, size_it);
    k = ceil(size_it / 2);
    [U,V, it, Error, GradientNorm] = LowRankAlgo(A, k, 10000, 1e-15, 1e-15, 'gradientnorm', 'eye', 'rand', 0);
    fprintf("Size %d x %d done\n", size_it, size_it);
    ThisExec = table(GradientNorm);
    CaseName = strcat(int2str(n_task), "_", int2str(size_it),'_grad_norm_analysis.csv');
    logPath = sprintf("./%s/%s", logsFolder, CaseName);
    writetable(ThisExec, logPath);
    size_it = size_it * 2;
    n_task = n_task + 1;
end