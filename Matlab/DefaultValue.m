function [v] = DefaultValue(arg)
switch arg
    case 'max_it'
        v = 1000;
    case 'err_eps'
        v = 1e-6;
    case 'grad_eps'
        v = 1e-6;
    case 'stop_c_type'
        v = 'approxerror';
    case 'init_t'
        v = 'eyedistinct';
    case 'init_w'
        v = 'zeros';
    case 'printsteps'
        v = 0;
    otherwise
        error('No such argument');
end
end
