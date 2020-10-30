% Executes the random restart version of the algorithm
% over a specific subset of the benchmark images 
% dataset and saves the statistics log.

datafolder = 'Frames_RandomRestart/abandonedBox';
logsFolder = 'Frames_RandomRestart_Out/abandonedBox';

if ~exist(logsFolder, 'dir')
       mkdir(logsFolder)
end

datapath = sprintf("%s/*.jpg", datafolder);
inputs = dir(datapath);

% Parameters of the algorithm
max_it = 1000;
err_eps = 1e-6;
grad_eps = 1e-6;
stop_c_type = 'default';
init_t = 'random';
init_w = 'default';
printsteps = 0;

% Number of random restarts
RRS = 10;

% Variable that expresses k as a percentage of the rank
k_perc = 0.3;

InputCount = length(inputs);

InputName = strings(InputCount, 1);
InputM = zeros(InputCount, 1);
InputN = zeros(InputCount, 1);
InputRk = zeros(InputCount, 1);
TargetRk = zeros(InputCount, 1);
Opt_delta = zeros(InputCount, 1);
Our_delta = zeros(InputCount, 1);
Delta_diff = zeros(InputCount, 1);
Our_Iter = zeros(InputCount, 1);


for fi = 1:InputCount
    
    FileName = inputs(fi).name;
    CaseName = erase(FileName, ".jpg");
    InputPath = sprintf("%s/%s", datafolder, FileName);
    
    % Converts the image to a matrix
    frame_img = imread(InputPath);
    gray_img = rgb2gray(frame_img);
    A = im2single(gray_img);
    
    [M,N] = size(A);
    Rk = rank(A);
    k = ceil(0.3 * Rk);
    [U,S,V] = svd(A);
    Ak_delta = norm(A- U(:,1:k)*S(1:k,1:k)*transpose(V(:,1:k)), 'fro') / norm(A, 'fro');
    
    fprintf("Task %d (%s): M=%d, N=%d, k=%d...\t", fi, CaseName, M, N, k);
    
    [U,V, it, Error, GradientNorm, RRErrors] = LowRankAlgoRR(RRS, A, k, max_it, err_eps, grad_eps, stop_c_type, init_t, init_w, printsteps);

    fprintf("...done!\n");
    
    delta = norm(A - U*V, 'fro')/norm(A, 'fro');
    
    InputName(fi) = CaseName;
    Opt_delta(fi) = Ak_delta;
    
    % Write the logs for this frame
    ThisExec = table(Error, GradientNorm);
    logPath = sprintf("%s/exec_%s.csv", logsFolder, CaseName);
    writetable(ThisExec, logPath);
    
    ThisRRs = table(RRErrors);
    logPathrr = sprintf("%s/rrs_%s.csv", logsFolder, CaseName);
    writetable(ThisRRs, logPathrr);
    
    clear('M', 'N', 'Rk', 'k', 'A', 'Ak_opt', 'Ak_delta', 'Time_opt');
end

% Write a summary log with the optimal errors
T = table(InputName, Opt_delta);
outPath = sprintf("%s/delta_opt.csv", execution_logs_folder);
writetable(T, outPath);