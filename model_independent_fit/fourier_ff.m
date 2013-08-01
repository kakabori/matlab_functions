function ret = fourierFF(params, qD)
%
%Use: ret = fourierFF(params, qD);
%
%INPUT
% params -> params(1) = F0, params(2) = F1 ...
% qD ->
%OUTPUT
% ret -> absolute value of form factor at qD values given by qD

if iscolumn(params) == 1
    params = transpose(parms);
end
if iscolumn(qD) == 1
    qD = transpose(qD);
end

ret = zeros(1,length(qD));
hmax = length(params)-1;
%First, calculate h = 0 term
%An applied constraint here is that rho(D/2)=rho_w, which determines F0 in terms of 
%a sum of Fh's: F0 = -2 * sum{Fh*(-1)^h} for h = 1 to hmax
%params(1) = 0;
%for h = 1 : hmax
%	params(1) = params(1) + params(h+1) * (-1)^h;
%end
%params(1) = -2 * params(1);
%ret = ret + params(1) * sin(0.5*qD) ./ (0.5*qD);

%Then, calculate h = 1 to hmax terms
for h = 1 : hmax
    ret = ret + params(h+1) * ( -2 * (-1)^h * sin(0.5*qD) ./ (0.5*qD) + sin(0.5*qD - pi*h) ./ (0.5*qD - pi*h) + sin(0.5*qD + pi*h) ./ (0.5*qD + pi*h));
end
ret = abs(ret);
