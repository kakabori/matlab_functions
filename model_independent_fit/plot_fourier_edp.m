function plotFourierEdp(D, params1, params2, params3)
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
if nargin > 2
    if iscolumn(params2) == 1
        params2 = transpose(params2);
    end
    if nargin > 3
        if iscolumn(params3) == 1
            params3 = transpose(params3);
        end
    end
end

zstep = 0.05;
z = -D/2:zstep:D/2;
zvec = z / D; %dimensionless z vector
Fh = params1 / D;
rho1 = fourierEdp(Fh, zvec);
if nargin > 2
    Fh = params2 / D;
    rho2 = fourierEdp(Fh, zvec);
    if nargin > 3
        Fh = params3 / D;
        rho3 = fourierEdp(Fh, zvec);
    end
end
figure; 
hold on;
hFit1 = plot(z, rho1);
if nargin > 2
    hFit2 = plot(z, rho2, '-g');
    if nargin > 3
        hFit3 = plot(z, rho3, '-c');
    end
end
hLegend = legend([hFit1, hFit2], 'Fit1', 'Fit2', 'location', 'NorthEast');