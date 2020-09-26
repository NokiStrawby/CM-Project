datafolder = 'Core_Random';

execution_logs_folder = 'logs';
logsFolder = sprintf('%s/%s', datafolder, execution_logs_folder);

if ~exist(logsFolder, 'dir')
       mkdir(logsFolder)
end

datapath = sprintf("%s/*.mat", datafolder);
inputs = dir(datapath);

max_it = 10000;
err_eps = 1e-6;
grad_eps = 1e-6;
stop_c_type = 'default';
init_t = 'default';
init_w = 'default';
printsteps = 0;

InputCount = length(inputs);

InputName = strings(InputCount, 1);
InputM = zeros(InputCount, 1);
InputN = zeros(InputCount, 1);
InputRk = zeros(InputCount, 1);
TargetRk = zeros(InputCount, 1);
Opt_delta = zeros(InputCount, 1);
Our_delta = zeros(InputCount, 1);
Delta_diff = zeros(InputCount, 1);
Svd_time = zeros(InputCount, 1);
Our_time = zeros(InputCount, 1);
Time_gain = zeros(InputCount, 1);
Our_Iter = zeros(InputCount, 1);


for fi = 1:InputCount
    FileName = inputs(fi).name;
    CaseName = FileName(1:end-4);
    InputPath = sprintf("%s/%s", datafolder, FileName);
    load(InputPath);
    
    fprintf("Task %d (%s): M=%d, N=%d, k=%d...\t", fi, CaseName, M, N, k);
    
    tic
    [U,V, it, Error, GradientNorm] = LowRankAlgo(A, k, max_it, err_eps, grad_eps, stop_c_type, init_t, init_w, printsteps);
    elapsed = toc;
    
    fprintf("...done!\n");
    
    delta = norm(A - U*V, 'fro')/norm(A, 'fro');
    % TODO: out step per convergenza
    
    InputName(fi) = CaseName;
    InputM(fi) = M;
    InputN(fi) = N;
    InputRk(fi) = Rk;
    TargetRk(fi) = k;
    Opt_delta(fi) = Ak_delta;
    Our_delta(fi) = delta;
    Delta_diff(fi) = delta-Ak_delta;
    Svd_time(fi) = Time_opt;
    Our_time(fi) = elapsed;
    Time_gain(fi) = Time_opt - elapsed;
    Our_Iter(fi) = it;
    
    ThisExec = table(Error, GradientNorm);
    logPath = sprintf("%s/exec_%s.csv", logsFolder, CaseName);
    writetable(ThisExec, logPath);
    
    clear('M', 'N', 'Rk', 'k', 'A', 'Ak_opt', 'Ak_delta', 'Time_opt');
end

T = table(InputName, InputM, InputN, InputRk, TargetRk, Opt_delta, Our_delta, Delta_diff, Svd_time, Our_time, Time_gain, Our_Iter);
T
outPath = sprintf("%s/summary.csv", datafolder);
writetable(T, outPath);