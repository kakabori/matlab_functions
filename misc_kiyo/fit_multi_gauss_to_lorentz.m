function [bestParams, chisquare] = fitMultiGauss2Lorentz(a0, gamma)
x = -1:0.01:1;
y = -1:0.01:1;
xdata = [];
ydata = [];

for k = 1 : length(x)
    xdata = [xdata x];
end
for k = 1 : length(y)
    temp = zeros(1, length(x)) + y(k);
    ydata = [ydata temp]; 
end

zdata = 1 ./ (xdata.^2 + ydata.^2 + (gamma/2)^2);
xydata = [xdata; ydata];

a0 = ones(1, 20); 
[bestParams, chisquare] = lsqcurvefit(@multiGauss, a0, xydata, zdata);

figure;
hold on;
plot(x, yLorentz(gamma, x, 0), 'b');
plot(x, yGauss(bestParams, x, 0), 'r');
plot(x, yLorentz(gamma, x, 0.2), 'b');
plot(x, yGauss(bestParams, x, 0.2), 'r');
plot(x, yLorentz(gamma, x, 0.4), 'b');
plot(x, yGauss(bestParams, x, 0.4), 'r');

end




function ret = multiGauss(a, data)
x = data(1, :);
y = data(2, :);
A = a(1:10);
s = a(11:20);
ret = zeros(1,length(x));

for k = 1 : 10    
    ret = ret + A(k)*exp(-(x.^2+y.^2)/(2*s(k)^2));
end
end

function ret = yGauss(a, x, y)
A = a(1:10);
s = a(11:20);
ret = zeros(1,length(x));

for k = 1 : 10 
    ret = ret + A(k)*exp(-(x.^2+y^2)/(2*s(k)^2));
end
end

function ret = yLorentz(gamma, x, y)
ret = 1 ./ (x.^2 + y^2 + (gamma/2)^2);
end