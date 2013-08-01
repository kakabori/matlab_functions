function [bestParams, chisquare] = fourierEdpFit(D, initial, phase, lb, ub, q, F, sigma)
%This program will perform nonlinear least square fit to 
%the input form factor by a model-independent method. 
%The bilayer elecron density profile is represented by 
%Fourier series.
%
%rho(z) = rho_w + (1/D) * sum{Fh * cos(2*pi*h*z/D)}
%F(q) = |sum{Fh * sinc(q*D/2-pi*h)}|,
%
%where sinc(x) = sin(x)/x and the sum runs from h = 0 to 
%hmax. Fh's are the fitting parameters. rho_w is the 
%electron density of bulk water. Input form factor is an 
%absolute value.
%
%EXAMPLE:
%[paramsfit,chisquare] = fourierEdpFit(D,initial,phase,[],[],q,F,sigma);
%[paramsfit,chisquare] = fourierEdpFit(D,initial,phase,lb,ub,q,F,sigma);
%
%INPUT
% D => D-spacing
% initial => row vector for initial parameter values. The 
%            first element corresponds to zeroth order
%            term, F0, and the last element to h_max-th 
%            order term, Fhmax.
% phase => row vector for the phases for each Fourier term.
%          -1 means the phase will be fixed to negative for 
%          the corresponding term and +1 positive. Assign 0
%          to fit the phase.
% lb => row vector for a set of lower bounds on parameters. 
%       Enter [] to use the default. The default is 
%       [-10 -10 -10 ...], the length equal to the length 
%       of initial vector
% ub => row vector for a set of upper bounds on parameters. 
%       Enter [] to use the default. The deualt is 
%       [10 10 10 ...], the length equal to the length of 
%       initial vector
% q => row vector for qz values
% F => row vector for F(qz) data points
% sigma => row vector for uncertainties on F(qz), sigma 
%          (not sigma^2)
%
%OUTPUT
% paramsfit -> output row vector for the best-fit parameters
% chisquare -> variable for chi squared value for the best fit
%
%DETAILS
% Let the order of Fourier series be hmax. This value is 
% determined by the following condition: hmax <= qmax*D/(2*pi), 
% which gives hmax = 8 for D = 62.8 Angstrom for form factor
% data that extends to 0.8 inverse Angstrom. To fit data with 
% 8th order Fourier series, which has 9 terms including zeroth 
% order, create four row vectors named initial, phase, lb, and 
% ub with 9 elements. For example, type, (white spaces are 
% entered for presentation purpose here)
%
% >> initial = [  1  -2  -2   2   2  -2  -2   1   1];
% >> phase   = [  0  -1  -1   1   1  -1  -1   0   0];
% >> lb      = [-10 -10 -10 -10 -10 -10 -10 -10 -10];
% >> ub      = [ 10  10  10  10  10  10  10  10  10];
%
% It is advised that a user always checks validity of input 
% vectors. Appropriate initial values, lower bounds, and upper 
% bounds should be chosen based on input phases. Make sure that 
% they respect the sign (positve/negative) of the phases and an 
% initial value is between its correponding parameter's lower 
% and upper bounds. See below for invalid inputs.
%
% Appropriate lower and upper bounds will be set based on the 
% input phases. For example, if the first element of initial 
% is 5, that of phase is -1, that of lb is 3, and that of ub 
% is 10, the lower bound on F0 will be set to -10 and the upper 
% bound on F0 will be set to -3. The inital value of F0 will 
% be set to -5. If abs(initial(1)) > 10 or abs(initial(1)) < 3, 
% the program will display an error message asking a user to 
% input an appropriate initial value for that particular element.

close all;

if D <= 0
    fprintf('D must be positive\n');
    return
end
for i = 1 : length(q)
    if q(i) < 0
        fprintf('q value must be zero or positive. Check q(%g)\n', i);
        return
    end
end
if isempty(lb)
    lb = zeros(1,length(initial));
    lb = lb - 10;
end
if isempty(ub)
    ub = zeros(1,length(initial));
    ub = ub + 10;
end

%Make sure all inputs have valid length
if length(initial) ~= length(phase)
    display('The number of elements in initial and phase vectors must be equal!\n');
    return
elseif length(initial) ~= length(lb)
    display('The number of elements in initial and lb vectors must be equal!\n');
    return
elseif length(initial) ~= length(ub)
    display('The number of elements in initial and ub vectors must be equal!\n');
    return
elseif length(q) ~= length(F)
    display('The number of elements in q and F vectors must be equal!\n');
    return
elseif length(q) ~= length(sigma)
    display('The number of elements in q and sigma vectors must be equal!\n');
    return
end

%Make sure all inputs are row vectors
if iscolumn(initial) == 1
    initial = transpose(initial);
end
if iscolumn(phase) == 1
    phase = transpose(phase);
end
if iscolumn(lb) == 1
    lb = transpose(lb);
end
if iscolumn(ub) == 1
    ub = transpose(ub);
end
if iscolumn(q) == 1
    q = transpose(q);
end
if iscolumn(F) == 1
    F = transpose(F);
end
if iscolumn(sigma) == 1
    sigma = transpose(sigma);
end

%xdata is the dimensionless version of q vector
%ydata is the form factor weighted by uncertainties, sigma
xdata = q * D;
ydata = F ./ sigma;

%Verify the input initial values, phases, and lower and upper bounds.
for a = 1 : length(initial)
    tempInitial = initial(a);
    tempPhase = phase(a);
    tempLb = lb(a);
    tempUb = ub(a);
    if tempLb > tempUb
        fprintf('lb(%g) > ub(%g). Lower bound cannot be greater than upper bound.\n', a, a);
        fpritnf('Check input lb and ub vectors.\n');
        return
    end
    if tempPhase == 1 %phase of ath term is positive, i.e., F_a >= 0
        initial(a) = abs(tempInitial); %make sure initial value is positive
        if tempLb < 0 && tempUb < 0
            fprintf('The upper bound must be positive for an order with a positive phase.\n');
            fprintf('Check consistency among phase(%g), lb(%g), and ub(%g)\n', a, a, a);
            return
        elseif tempLb < 0 %set correct lower bound. upper bound is already correct.
            lb(a) = 0;
        end
        if initial(a) > ub(a) || initial(a) < lb(a)
            fprintf('The initial value must be between the lower and upper bounds\n');
            fprintf('Check consistency among initial(%g), lb(%g), and ub(%g)\n', a, a, a);
            return
        end
    elseif tempPhase == -1 %phase of ath term is negative
        initial(a) = -abs(tempInitial); %make sure initial value is negative
        if tempLb > 0 && tempUb > 0 
            fprintf('The lower bound must be negative for an order with a negative phase.\n');
            fprintf('Check consistency among phase(%g), lb(%g), and ub(%g)\n', a, a, a);
            return
        elseif tempUb > 0 %set correct upper bound. lower bound is already correct.
            ub(a) = 0;
        end
        if initial(a) > ub(a) || initial(a) < lb(a)
            fprintf('The initial value must be between the lower and upper bounds\n');
            fprintf('Check consistency among initial(%g), lb(%g), and ub(%g)\n', a, a, a);
            return
        end
    elseif tempPhase == 0 %phase of ath term is a free variable
        if tempInitial > tempUb || tempInitial < tempLb
            fprintf('The initial value must be between the lower and upper bounds\n');
            fprintf('Check consistency among initial(%g), lb(%g), and ub(%g)\n', a, a, a);
            return
        end
    else
        fprintf('The phases must be 0, -1, or 1. Check phase(%g).\n', a);
        return
    end
end

myfunc = @(initial, xdata)(fourierFF(initial, xdata) ./ sigma);
[bestParams, chisquare] = lsqcurvefit(myfunc, initial, xdata, ydata, lb, ub);
%[paramsfit, chisquare] = lsqcurvefit(@scalingFactor, params, xdata, ydata, lb, ub, options);

%Due to a constraint that rho(D/2)=rho_w, F0 is determined in terms of 
%the sum of Fh's: F0 = -2 * sum{Fh*(-1)^h} for h = 1 to hmax
bestParams(1) = 0;
hmax = length(bestParams)-1;
for h = 1 : hmax
	bestParams(1) = bestParams(1) + bestParams(h+1) * (-1)^h;
end
bestParams(1) = -2 * bestParams(1);

%Now, calculate and show the fitted form factor with the data
qz = 0:0.002:1.0;
qzvec = qz * D;
yfit = fourierFF(bestParams, qzvec);
figure; 
hold on;
hE = errorbar(q, F, sigma, 'xr');
hFit = plot(qz, yfit);
errorbar_tick(hE, 1000);

%Calculate electron density profile
zstep = 0.05;
z = -D/2:zstep:D/2;
zvec = z / D; %dimensionless z vector
Fh = bestParams / D;
rho = fourierEdp(Fh, zvec);
figure; plot(z, rho);
