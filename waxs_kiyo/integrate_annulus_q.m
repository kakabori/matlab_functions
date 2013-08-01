function result = integrate_annulus_q(imag, q_range, phi_range, bin_size, N_point, q_delta)
%integrate_annulus_q
%This function is integrate_annulus for an image that is already transformed 
%from CCD to q-space. 'imag' must be an output from transform_ccd2q
%function, which follows a specific format. See transform_ccd2q.m for the
%details of the format.
%
%   imag ->
%   q_range->
%   phi_range ->
%   bin_size -> size of angular bin in degree.
%   N_point -> number of interpolation point in phi within each bin.
%   q_delta -> resolution of interpolation in q in inverse Angstrom.

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

%Extract qr and qz labels, and the intensity matrix according to the format
%specified by KA. See transform_ccd2q.m regarding the format.
[m n] = size(imag);
qr = imag(1,2:n);
qz = imag(2:m,1);
Int = imag(2:m,2:n);
delta_qz = qz(2)-qz(1);
delta_qr = qr(2)-qr(1);

%Set q_delta to the smaller of delta_qz and delta_qr
if (nargin < 6)
    if (delta_qz <= delta_qr)
        q_delta = delta_qz;
    else
        q_delta = delta_qr;
    end;
end;

%beam_col = find(qr == 0);These probably are not necessary.
%beam_row = find(qz == 0);

phi_delta = bin_size/N_point;
phi_min = phi_range(1)-phi_delta*(N_point-1)/2;
phi_max = phi_range(2)+phi_delta*(N_point-1)/2;

%Set up the matrices whose elements are the values of q and phi for interpolation. 
q = q_range(1):q_delta:q_range(2);
phi = deg2rad(phi_min):deg2rad(phi_delta):deg2rad(phi_max);
len_q = length(q);
len_phi = length(phi);
q2 = repmat(q,[len_phi 1]);
phi2 = repmat(phi',[1 len_q]);

%Construct the matrices that will be used as the grid in interpolation.
len_qr=length(qr);
len_qz=length(qz);
qr2=repmat(qr,[len_qz,1]);
qz2=repmat(qz,[1,len_qr]);

%Calculate the values of (qr,qz) at which intensity will be interpolated.
XI = q2 .* cos(phi2);
YI = q2 .* sin(phi2);

%Do the interpolation.
I = interp2(qr2,qz2,Int,XI,YI,'spline');

%Find the points that are outside of the input image.
A = find(XI > max(qr));
B = find(XI < min(qr));
C = find(YI > max(qz));
D = find(YI < min(qz));
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