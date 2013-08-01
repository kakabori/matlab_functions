function ret = fourierEdp(params, zvec)
%This function returns a row vector, rho.
% params -> Fh / D
% zvec -> z / D

%rho_w = 0.33333;

if iscolumn(params) == 1
    params = transpose(params);
end
if iscolumn(zvec) == 1
    zvec = transpose(zvec);
end

ret = zeros(1, length(zvec));
hmax = length(params)-1;
for h = 1 : hmax
    ret = ret + params(h+1) * ((-1)^(h+1) + cos(2*pi*h*zvec));
end
ret = 2 * ret;
%ret = ret + rho_w;
