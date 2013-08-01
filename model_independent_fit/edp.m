function rho = edp(params, zvec)
%This function returns a row vector, rho.
%   Detailed explanation goes here
D = 62.9;
if isrow(zvec) == 1
    zvec = transpose(zvec);
end
matrix = zeros(length(zvec),length(params));
for k = 1 : length(params)
    matrix(:,k) = params(k) * ((-1)^(k+1) + cos(2*pi*k*zvec/D));
end
rho = sum(matrix,2);
rho = transpose(rho);
%figure; plot(zvec,rho);
end
