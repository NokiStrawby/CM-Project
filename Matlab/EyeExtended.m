% Returns a matrix of size m x n shaped as (suppose m > n)
%{
    |1 0 0 1 0|
    |0 1 0 0 1|
    |0 0 1 0 0|
%}
% While if distinct = true we generate a similar shape but with
% non-duplicated columns, as:
%{
    |1 0 0 2 0 0 3|
    |0 1 0 0 2 0 0|
    |0 0 1 0 0 2 0|
%}

function [X] = EyeExtended(m,n, distinct)
arguments % Parsing input arguments
    m (1,1) {mustBeNumeric}
    n (1,1) {mustBeNumeric}
    distinct = false
end

r = min(m,n);
c = max(m,n);

mr = r;
t = ceil(c/r);
mc = r * t;
M = zeros(mr, mc);

for i = 1:t
    cs = (i-1)*mr +1;
    ce = i*mr;
    toSet = eye(mr);
    if distinct
        toSet = i*toSet;
    end
    M(:,cs:ce) = toSet; 
end

if (m > n)
    M = M';
end

X = M(1:m, 1:n);
end

