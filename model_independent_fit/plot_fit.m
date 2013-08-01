function plotFit(params, xdata, ydata)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if iscolumn(xdata) == 1
    xdata = transpose(xdata);
end
if iscolumn(ydata) == 1
    ydata = transpose(ydata);
end
yfit = scalingFactor(params, xdata);
figure; plot(xdata, ydata, '+r', xdata, yfit, 'b');

temp = yfit - ydata;
temp = temp .* temp;
rms = sum(temp);
fprintf('rms = %g\n', rms);

end

