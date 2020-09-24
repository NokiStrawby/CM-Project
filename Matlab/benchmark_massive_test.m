listDir = dir("cdnet14/dataset/*");
for sub_dir_idx = 3:size(listDir,1)
    subDir_name = listDir(sub_dir_idx).name;
    subDir_path = strcat("cdnet14/dataset/", subDir_name); 
    subDir_frames = dir(strcat(subDir_path, "/input/*.jpg"));
    sz = size(subDir_frames,1);
    frame_step = ceil(sz/30);
    dir_res = cell(sz*20,10);
    dir_res_it = 0;
    dir_tic = tic;
    k_perc = 0.3;

    strcat("Directory ", subDir_name)
    for frame_idx = 1:frame_step:size(subDir_frames,1)
        out_path = strcat(subDir_path, "/out");
        [status, msg, msgID] = mkdir(out_path);
        frame_name = subDir_frames(frame_idx).name;
        strcat("Frame ", frame_name)
        frame_path = strcat(subDir_path + "/input/" + frame_name);
        frame_img = imread(frame_path);
        gray_img = rgb2gray(frame_img);
        img_matrix = im2single(gray_img);
        img_matr_rank = rank(img_matrix);
        [U,S,V] = svd(img_matrix);
        [m,n] = size(img_matrix);
        k = ceil(k_perc * img_matr_rank);
        matr_data = {erase(frame_name, ".jpg"), m, n, img_matr_rank};
        opt_approx_err = norm(img_matrix - U(:,1:k)*S(1:k,1:k)*transpose(V(:,1:k))) / norm(img_matrix);
        [err, diff, el_time, it] = ApproxFrames(img_matrix, k, opt_approx_err);
        new_row = [matr_data, k, err, opt_approx_err,  diff, el_time, it];
        dir_res_it = dir_res_it + 1;
        dir_res(dir_res_it, 1:10) = new_row;
        strcat("Frame ", frame_name, " k ", string(k), " completato ", string(el_time))
    end
    dir_res = dir_res(1:dir_res_it, 1:10);
    T = cell2table(dir_res, 'VariableNames', {'frame_name', 'm', 'n', 'rank', 'k', 'err', 'opt', 'diff', 'el_time', 'it'});
    writetable(T, strcat(subDir_path, "/out/frames_analysis_", subDir_name, ".csv"));
    strcat("Directory ", subDir_name, " completata ", string(toc(dir_tic)))
end

function [err, diff, el_time, it] = ApproxFrames(A, k, opt_err)
    frame_k_tic = tic;
    [U,V,it,~,~] = LowRankAlgo(A, k, 1000, 1e-5, 1e-15, 'approxerror', 'eyedistinct', 'zeros', 0);
    el_time = toc(frame_k_tic);
    err = norm(A - U*V) / norm(A);
    diff = err - opt_err;
end