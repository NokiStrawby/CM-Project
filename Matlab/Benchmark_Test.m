% Executes the algorithm over (a certain subset of) the benchmark 
% frames dataset and saves the execution stats

data_Dir = "cdnet14/dataset/*";
out_Dir = "cdnet14/out";
sub_Dirs = dir(data_Dir);

if ~exist(out_Dir, 'dir')
       mkdir(out_Dir)
end

% Variable that determines on how many frames the algorithm
% has to be executed for each sub-directory
n_frames = 100;

% Variable that expresses k as a percentage of the rank
k_perc = 0.3;

for sub_dir_idx = 3:size(sub_Dirs,1)
    
    subDir_name = listDir(sub_dir_idx).name;
    subDir_path = strcat("cdnet14/dataset/", subDir_name); 
    subDir_frames = dir(strcat(subDir_path, "/input/*.jpg"));
    
    sz = size(subDir_frames,1);
    % Set to 1 for an execution on all frames
    frame_step = ceil(sz/n_frames);
    
    res_tab = cell(sz,10);
    res_it = 0;

    for frame_idx = 1:frame_step:sz
        
        frame_name = subDir_frames(frame_idx).name;
        frame_path = strcat(subDir_path + "/input/" + frame_name);
        
        % Converts the image to a matrix
        frame_img = imread(frame_path);
        gray_img = rgb2gray(frame_img);
        img_matrix = im2single(gray_img);
        
        img_matr_rank = rank(img_matrix);
        k = ceil(k_perc * img_matr_rank);
        
        [U,S,V] = svd(img_matrix);
        [m,n] = size(img_matrix);
        
        matr_data = {erase(frame_name, ".jpg"), m, n, img_matr_rank};
        trunc_svd_approx = U(:,1:k)*S(1:k,1:k)*transpose(V(:,1:k));
        
        opt_approx_err = norm(img_matrix - trunc_svd_approx, 'fro') / norm(img_matrix, 'fro');
        
        [err, diff, el_time, it] = ApproxFrame(img_matrix, k, opt_approx_err);
        
        res_it = res_it + 1;
        new_row = [matr_data, k, err, opt_approx_err,  diff, el_time, it];
        res_tab(res_it, 1:10) = new_row;
        
        strcat("Frame ", frame_name, " completato ", string(el_time))
    end
    
    res_tab = res_tab(1:res_it, 1:10);
    
    T = cell2table(res_tab, 'VariableNames', {'frame_name', 'm', 'n', 'rank', 'k', 'err', 'opt', 'diff', 'el_time', 'it'});
    writetable(T, strcat(out_Dir, "/frames_analysis_", subDir_name, ".csv"));
    strcat("Directory ", subDir_name, " completata")
end

% Approximates a frame matrix with the Low Rank Approximation Algorithm
% and returns the approximation error, the difference to the optimal error,
% the elapsed time and the number of iterations
function [err, diff, el_time, it] = ApproxFrame(A, k, opt_err)

    frame_k_tic = tic;
    [U,V,it,~,~] = LowRankAlgo(A, k, 1000, 1e-6, 1e-15, 'approxerror', 'eye', 'zeros', 0);
    el_time = toc(frame_k_tic);
    
    err = norm(A - U*V, 'fro') / norm(A, 'fro');
    diff = err - opt_err;
end