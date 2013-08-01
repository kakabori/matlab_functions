function ret = calcChi2(gamma, params, xdata, ydata)
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

theory = yGauss(params, xdata, y);
diff = (theory - ydata) ./ sigma;
diff = diff .* diff;
ret = sum(diff);
end

function ret = yGauss(a, x, y)
global kmax;
A = a(1:kmax);
s = a((kmax+1):(2*kmax));
z0 = a(2*kmax+1);
m = a(2*kmax+2);
m = -0.1;
ret = zeros(1,length(x));
%ret = ret + z0;
ret = ret + z0 + m*y;

for k = 1 : kmax
    ret = ret + A(k)*exp(-(x.^2+y^2)/(2*s(k)^2));
end
end

function ret = yLorentz(gamma, x, y)
ret = (gamma / pi) ./ (x.^2 + y^2 + gamma^2);
end