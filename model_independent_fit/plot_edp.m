function plotEdp(params, z, data)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
rho = edp(params, z);
figure; plot(z, data, '+r', z, rho, 'b');
end

