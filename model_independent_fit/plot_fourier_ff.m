function plotFourierFF(D, xdata, ydata, sigma, params1, params2, params3)
%
%Use: ret = fourierFF(params, qD);
%
%INPUT
% params ->
% qD ->
%OUTPUT
% ret -> absolute value of form factor at qD values given by qD

if iscolumn(params1) == 1
    params1 = transpose(params1);
end
if nargin > 5
    if iscolumn(params2) == 1
        params2 = transpose(params2);
    end
    if nargin > 6
        if iscolumn(params3) == 1
            params3 = transpose(params3);
        end
    end
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

qz = 0:0.002:1.0;
qzvec = qz * D;
yfit1 = fourierFF(params1, qzvec);
if nargin > 5
    yfit2 = fourierFF(params2, qzvec);
    if nargin > 6
        yfit3 = fourierFF(params3, qzvec);
    end
end
figure; 
hold on;
hE = errorbar(xdata, ydata, sigma, 'xr');
hFit1 = plot(qz, yfit1);
if nargin > 5
    hFit2 = plot(qz, yfit2, '-g');
    if nargin > 6
        hFit3 = plot(qz, yfit3, '-c');
    end
end
errorbar_tick(hE, 1000);
hLegend = legend([hE, hFit1, hFit2], 'Data', 'Fit1', 'Fit2', 'location', 'NorthEast');