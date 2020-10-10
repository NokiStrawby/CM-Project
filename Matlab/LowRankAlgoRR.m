% Executes the algorithm with a random restart approach, having a number of
% random restarts equal to rrs

% The parameter init_t must be set to use a random initialization logic
% (that is, one among random, randfull and randeye)

function [U,V, it, errs, gradnorms, rrerrors] = LowRankAlgoRR(rrs, A, k, max_it, err_eps, grad_eps, stop_c_type, init_t, init_w, printsteps)
arguments
    rrs (1,1) {mustBeNumeric}
    A (:,:) {mustBeNumeric}
    k (1,1) {mustBeNumeric}
    max_it (1,1) {mustBeNumeric} = DefaultValue('max_it')
    err_eps (1,1) {mustBeNumeric} = DefaultValue('err_eps')
    grad_eps (1,1) {mustBeNumeric} = DefaultValue('grad_eps')
    stop_c_type = DefaultValue('stop_c_type')
    init_t = 'default'
    init_w = DefaultValue('init_w')
    printsteps = DefaultValue('printsteps')
end

if strcmp(init_t, 'default')
    init_t = 'random';
end

randomInits = {'random', 'randfull', 'randeye'};

if ~any(strcmp(randomInits, init_t))
    % Check that the initialization mus be random
    error('Initialization logic cant be constant')
end

best_err = realmax();
store_errors = zeros(rrs,1);

for i = 1:rrs
    % Execute the algorithm
    [U_i,V_i, it_i, errs_i, gns_i] = LowRankAlgo(A, k, max_it, err_eps, grad_eps, stop_c_type, init_t, init_w, printsteps);
    
    if errs_i(end) < best_err
        best_U = U_i;
        best_V = V_i;
        best_it = it_i;
        best_errs = errs_i;
        best_gns = gns_i;
        best_err = errs_i(end);
    end
    
    store_errors(i) = errs_i(end);
end

U = best_U;
V = best_V;
it = best_it;
errs = best_errs;
gradnorms = best_gns;
rrerrors = store_errors;

end

