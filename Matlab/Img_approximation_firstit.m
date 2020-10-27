im_path = "Img_Approx/cat.png";
out_path = "Img_Approx/out/cat_firstit";

im_in = imread(im_path);
gray_im = rgb2gray(im_in);
im_matr = im2single(gray_im);
rk = rank(im_matr);
k = ceil(0.025 * rk);

tic
[U_svd, S, V_svd] = svd(im_matr);
fprintf("Elapsed time SVD %f\n", toc);

for it=1:10
   tic
   [U,V,~,~,~] = LowRankAlgo(im_matr, k, it, 1e-6, 1e-15, 'approxerror', 'eyeextended', 'zeros', 0);
   fprintf("Numero iterazioni %d, elapsed time %f\n", it, toc);
   alg_appr = U*V;
   alg_err = norm(im_matr - alg_appr, 'fro') / norm(im_matr, 'fro');
   out_alg_name = sprintf("%s/res_lralg_%d.png",out_path,it);
   imwrite(alg_appr, out_alg_name);   
end