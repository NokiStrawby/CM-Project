% Executes the algorithm multiple times over a specified matrix 
% derived from an image, for a specific value of k.
% For each execution, stop_it is set for it to stop after a specific
% iteration. The resulting approximation matrix of each execution
% is saved as an image.

im_path = "Img_Approx/cat.png";
out_path = "Img_Approx/out/cat_firstit";

% Variable that expresses k as a percentage of the rank
k_perc = 0.5;

% Variable that expresses the number of times the algorithm has to be 
% executed. Each execution is stopped after an increasing number of
% iterations.
num_it = 10;

% Converts the image to a matrix
im_in = imread(im_path);
gray_im = rgb2gray(im_in);
im_matr = im2single(gray_im);
rk = rank(im_matr);

k = ceil(k_perc * rk);

tic
[U_svd, S, V_svd] = svd(im_matr);
fprintf("Elapsed time SVD %f\n", toc);


for it=1:num_it
   tic
   [U,V,~,~,~] = LowRankAlgo(im_matr, k, it, 1e-6, 1e-15, 'approxerror', 'eyeextended', 'zeros', 0);
   fprintf("Iterations %d, elapsed time %f\n", it, toc);
   
   alg_appr = U*V;
   alg_err = norm(im_matr - alg_appr, 'fro') / norm(im_matr, 'fro');
   out_alg_name = sprintf("%s/res_lralg_%d.png",out_path,it);
   imwrite(alg_appr, out_alg_name);   
end