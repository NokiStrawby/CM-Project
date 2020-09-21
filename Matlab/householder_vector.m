function [v, s] = householder_vector(x) 
if x(1) == 0
    s = 1;
else
    s = -sign(x(1));
end
s = s * norm(x);

v = x;
v(1) = v(1) - s;
v = v / norm(v);