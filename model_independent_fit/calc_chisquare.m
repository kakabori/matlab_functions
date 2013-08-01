function ret = calcChisquare(D, params, xdata, ydata, sigma)
%
%Use: ret = fourierFF(params, qD);
%
%INPUT
% params ->
% qD ->
%OUTPUT
% ret -> absolute value of form factor at qD values given by qD

if iscolumn(params) == 1
    params = transpose(parms);
end
if iscolumn(xdata) == 1
    xdata = transpose(xdata);
end
if iscolumn(ydata) == 1
    ydata = transpose(ydata);
end
if iscolumn(sigma) == 1
    sigma = transpose(sigma);
end

qD = xdata * D;

theory = fourierFF(params, qD);
diff = (theory - ydata) ./ sigma;
diff = diff .* diff;
ret = sum(diff);