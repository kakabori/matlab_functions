function [bestParams, chisquare, gset] = fitGauss2Lorentz(gamma1, gamma2, a0)
%[bestParams, chisquare, gset] = fitGauss2Lorentz(gamma1, gamma2, a0)
%This function performs nonlinear least square fit to 2D Lorentzian with seven Gaussians and one
%constant, totaling 15 free parameters. It performs the fit at 5000 different angles. The range of 
%angles is given by gamma1 and gamma2. a0 is a 15-element row vector that contains initial values.
%a(1,1) to a(1,7) are initial values for Gaussain amplitudes and a(1,8) to a(1,14) are initial
%values for Gaussian sigma. a(1,15) is an initial value for a constant offset.
%

tic
global kmax;
kmax = 7;

if nargin == 2
    a0 = ones(1,2*kmax+1);
    %a0(1,1) = 50;
    %a0(1,2) = 50;
    a0(1,2*kmax+1) = 0.001;
end
lb = zeros(1, 2*kmax+1);
lb(1,2*kmax+1) = -10;
ub = 7000 * ones(1, 2*kmax+1);
%options = optimset('MaxFunEvals', 10000, 'TolFun', 10^-8, 'TolX', 10^-8);
options = optimset('Display', 'off', 'MaxFunEvals', 10000);

g1 = log10(gamma1);
g2 = log10(gamma2);
step = (g2 - g1) / 5000;
if step == 0
    step = 1;
end
g12 = g1:step:g2;
gset = 10.^g12;
bestParams = [];
chisquare = [];
for g = gset
    x = (-50*g):(g/10):(50*g);
    y = zeros(1,length(x));
    xydata = [x; y];
    zdata = (g / pi) ./ (x.^2 + y.^2 + g^2);
    [best, chi] = lsqcurvefit(@multiGauss, a0, xydata, zdata, lb, ub, options);
    bestParams = [bestParams; best];
    chisquare = [chisquare; chi];
    a0 = best;
end

figure; hold on;
for k = 1 : kmax
    plot(gset, bestParams(:,k));
end
figure; hold on;
for k = (kmax+1) : 2*kmax
    plot(gset, bestParams(:,k));
end
figure; 
plot(gset, bestParams(:,(2*kmax+1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%probably unnecessary
% figure;
% plot(bestParams(:,(2*kmax+2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% hold on;
% plot(x, yLorentz(gamma1, x, 0*gamma1), 'b');
% plot(x, yGauss(bestParams(1,:), x, 0*gamma1), 'r');
% plot(x, yLorentz(gamma1, x, 5*gamma1), 'b');
% plot(x, yGauss(bestParams(1,:), x, 5*gamma1), 'r');
% plot(x, yLorentz(gamma1, x, 10*gamma1), 'b');
% plot(x, yGauss(bestParams(1,:), x, 10*gamma1), 'r');

gset = gset';
toc
end




function ret = multiGauss(a, data)
global kmax;
x = data(1, :);
y = data(2, :);
A = a(1:kmax);
s = a((kmax+1):(2*kmax));
z0 = a(2*kmax+1);
ret = zeros(1,length(x));
ret = ret + z0;

for k = 1 : kmax    
    ret = ret + A(k)*exp(-(x.^2+y.^2)/(2*s(k)^2));
end
end

function ret = yGauss(a, x, y)
global kmax;
A = a(1:kmax);
s = a((kmax+1):(2*kmax));
z0 = a(2*kmax+1);
ret = zeros(1,length(x));
ret = ret + z0;

for k = 1 : kmax
    ret = ret + A(k)*exp(-(x.^2+y^2)/(2*s(k)^2));
end
end

function ret = yLorentz(gamma, x, y)
ret = (gamma / pi) ./ (x.^2 + y^2 + gamma^2);
end