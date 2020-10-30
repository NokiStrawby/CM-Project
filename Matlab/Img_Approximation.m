im_path = "Img_Approx/cat.png";
out_path = "Img_Approx/out/cat";

im_in = imread(im_path);
gray_im = rgb2gray(im_in);
im_matr = im2single(gray_im);
rk = rank(im_matr);
[U_svd,S,V_svd] = svd(im_matr);
ks = [0.6*rk; 0.2*rk; 0.05*rk; 0.025*rk; 0.01*rk];


for k_i = 1:5
   k = ceil(ks(k_i,1));
   fprintf("k = %d\n", k);
   opt_appr = U_svd(:,1:k)*S(1:k,1:k)*transpose(V_svd(:,1:k));
   opt_err = norm(im_matr - opt_appr, 'fro') / norm(im_matr, 'fro');
   [U,V,~,~,~] = LowRankAlgo(im_matr, k, 2000, 1e-6, 1e-15, 'approxerror', 'eye', 'zeros', 0);
   alg_appr = U*V;
   alg_err = norm(im_matr - alg_appr, 'fro') / norm(im_matr, 'fro');
   alg_err - opt_err
   out_svd_name = sprintf("%s/res_svd_%d.png",out_path,k);
   out_alg_name = sprintf("%s/res_lralg_%d.png",out_path,k);
   imwrite(opt_appr, out_svd_name);
   imwrite(alg_appr, out_alg_name);   
   fprintf('%d finito\n ', k);
end






