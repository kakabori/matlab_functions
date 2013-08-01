function [ chisquare parameters ] = edpfit_recursive(num, xdata, ydata)
%UNTITLED3 Summary of this function goes here
% num -> number of parameters
% xdata -> qz vector
% ydata -> scaling factor vector
%% Initial setup
% Define global variables used for fitting
global largeNum;
largeNum = 10^10;
% Set maximum number of recursion (iteration)
maxRecur = 1000;
% Make sure that xdata and ydata are both row vectors
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
% Create a holder for chisquare values. This is a column vector
chisquare = [];
% Create a holder for fitted parameter values. This is a matrix. Each row contains fitted parameter
% values.
parameters = [];
% A variable to holder the best chisquare
best = 0;
%% Start recursive fitting procedure
for k = 1 : maxRecur
    % First, randomly determine initial parameter values
    params = random('Uniform', -10, 10, 1, num);
    [temppara, tempchi] = lsqcurvefit(@scalingFactor, params, xdata, ydata, [], [], options);
    if tempchi <= largeNum^2
        chisquare = [chisquare; tempchi];
        parameters = [parameters; temppara];
    end
end

end

