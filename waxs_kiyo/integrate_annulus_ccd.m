function result = integrate_annulus_ccd(imag, q_range, phi_range, bin_size, N_point, q_delta)
%result = integrate_annulus_ccd(imag, q_range, phi_range, bin_size, N_point, q_delta)
%
%Function to perform annular integral on x-ray scattering image
%
%INPUT
% imag: N*M x-ray scattering image
% q_range: Integrate annulus for q_range(1) < q < q_range(2)
%          Units are reciprocal Angstroms
%          Default values are  1.0 A^-1 to 1.6 A^-1
% phi_range: Calculate for phi_range(1) < phi < phi_range(2)
%            Angles are in degrees and measured
%            counter-clockwise from (+) horizontal (right)
%            Defaults are 0 degrees to 90 degrees
% bin_size: 
%           Default value is 1 degree
% N_point: number of interpolation point in phi within each bin
%          Default value is 1
% q_delta: resolution of interpolation in q in inverse Angstrom
%          Default value is 0.0025
%
%OUTPUT
% result: 2D matrix
%
%DETAILS
%

% Get Global variables
global X_cen Y_cen; global X_Lambda Spec_to_Phos;

% Set q_range 
if (nargin<2)
    q_range = [1.0 1.6]; % Default q-range is 1.0 A^-1 < q < 1.6 A^-1.
end

% Set phi_range
if (nargin<3)
    phi_range(1) = 0;  % Start at 0 degrees (horizontal axis going right)
    phi_range(2) = 90; % Finish at 90 degrees (vertical axis going up).
end

% Set size of angular bin in degrees.
if (nargin<4)
    bin_size = 1;
end

% Set number of interpolation points in phi.
if (nargin<5)
    N_point = 1;
end

% Set resolution of interpolation in q in inverse Angstrom.
if (nargin<6)
    q_delta = 0.0025;
end

phi_delta = bin_size/N_point;
phi_min = phi_range(1)-phi_delta*(N_point-1)/2;
phi_max = phi_range(2)+phi_delta*(N_point-1)/2;

q = q_range(1):q_delta:q_range(2);
phi = deg2rad(phi_min):deg2rad(phi_delta):deg2rad(phi_max);
len_q = length(q);
len_phi = length(phi);

q2 = repmat(q,[len_phi 1]);
phi2 = repmat(phi',[1 len_q]);

theta = asin(q2*X_Lambda/4/pi);
R = Spec_to_Phos*tan(2*theta);
sin_phi = sin(phi2);
cos_phi = cos(phi2);
y = R.*sin_phi;
x = R.*cos_phi;

X = Y_cen + x;
Y = X_cen - y;

I = interp2(imag,X,Y,'spline');
%figure; colormap gray;
%imagesc(q,rad2deg(phi),I,[0 6000]);

A = find(X > 1024);
B = find(X < 1);
C = find(Y > 1024);
D = find(Y < 1);
I(A) = 0;
I(B) = 0;
I(C) = 0;
I(D) = 0;

% Sum intensity along q's for each phi and create a column vector, I, as a
% function of phi
I_phi = sum(I,2)/length(q);

% Calculate number of bins.
N_bin = length(I_phi)/N_point;

% Sum intensity within each bin and create a reduced column vector, I_phi2.
I_phi2 = zeros(1,N_bin);
for j = 1:N_bin
    I_phi2(j) = sum(I_phi(((j-1)*N_point+1):j*N_point));
end
I_phi2 = I_phi2/N_point;

% Transpose I_phi2 vector to make it a column vector
I_phi2 = I_phi2';

% Create the label for phi angles.
phi_bin = phi_range(1):bin_size:phi_range(2);
phi_bin = phi_bin';

% Create the output
result = [phi_bin I_phi2];

end
