function [paramsfit, chisquare] = edpfit(params, xdata, ydata)
% params->
% xdata->
% ydata->

close all;
global largeNum;
largeNum = 10^10;

if iscolumn(xdata) == 1
    xdata = transpose(xdata);
end
if iscolumn(ydata) == 1
    ydata = transpose(ydata);
end
% Add another element to xdata and ydata. The additional element will be used to give huge chisquare
% to an unreasonable fit. 
xdata = [xdata 100];
ydata = [ydata 0];
options = optimset('MaxFunEvals', 1000000, 'MaxIter', 2000, 'Display', 'off');

[paramsfit, chisquare] = lsqcurvefit(@scalingFactor, params, xdata, ydata, [], [], options);
yfit = scalingFactor(paramsfit, xdata);
figure; plot(xdata(1:end-1), ydata(1:end-1), '+r', xdata(1:end-1), yfit(1:end-1), 'b');

%% Calculate electron density profile
D = 62.9;
zstep = 0.5;
zvec = [-D/2:zstep:D/2];
zvec = transpose(zvec);
matrix = zeros(length(zvec),length(paramsfit));
for k = 1 : length(paramsfit)
    matrix(:,k) = paramsfit(k) * ((-1)^(k+1) + cos(2*pi*k*zvec/D));
end
rho = sum(matrix,2);
% check boundary condition
temp = 0;
for k = 1 : length(paramsfit)
    temp = temp + paramsfit(k) * (-1)^(k-1);
end


zvec = transpose(zvec);
rho = transpose(rho);
figure; plot(zvec,rho);
end

