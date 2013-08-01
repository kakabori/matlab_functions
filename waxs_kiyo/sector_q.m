function result = sector_q(imag, q_range, phi_range, q_delta, phi_delta)
%
%This function will return intensity as a function of q along a constant value
%of phi on a q-corrected image. The output will contain two column vectors, 
%q and I. 
%The program will first transform the input data in 
%Cartensian coordinates into polar coordinates by interpolating within a region 
%specified by q_range and phi_range. Data points will be created at
%every q_delta and phi_delta. 
%Then, the intensity will be averaged over a range specified by phi_range. 
%Data points will finally be created at every q_delta for a range specified 
%by q_range. 
%If no averaging is desired, phi_range(1) should be set equal to phi_range(2).
%
%Caveat: Since the input image will be interpolated, a user should pay
%attention to the resolution of the image and avoid
%interpolating many points between two real data points. The resolution in
%phi is sensitive to the position on the image; very close to the origin,
%phi resolution is very coarse.
%   
%   Inputs  
%       imag -> N*M x-ray scattering image in Cartesian coordinates in q-space.
%       The input image must be an output of transform_ccd2q function so
%       that it follows the specific format defined by KA. See
%       transform_ccd2q.m for the details of the format.
%
%       q_range -> create plot for q_range(1) < q < q_range(2).
%       Units are inverse Angstroms.
%       
%       phi_range -> will integrate intesity for phi_range(1) < phi < phi_range(2).
%       Angles are in degrees and measured counter-clockwise 
%       from the positive equator.
%                    
%       q_delta -> Resolution in q for the plot. Units are inverse Angstroms.
%       
%       phi_delta -> Finess of interpolation in phi. Units are degrees.
%       After interpolation,
%       the program will average the intensity over phi_range. phi_delta
%       should not be much smaller than the resolution of the input image
%       itself.
% 
%   Output 
%       result = [q, I], where both q and I are column vectors. Use these 
%       variables to plot intensity as a function of q in MATLAB.
%
%   KA 3/14/2012. 

% Set phi_range
if (nargin<3)
    phi_range(1) = 0;  % Start at 0 degrees (horizontal axis going right)
    phi_range(2) = 90; % Finish at 90 degrees (vertical axis going up).
end

% Set resolution in q
if (nargin<4)
    q_delta = 0.1;
end

% Set the step size of interpolation in phi.
if (nargin<5)
    phi_delta = 0.1;
end

%Extract qr and qz labels, and the intensity matrix according to the format
%specified by KA. See transform_ccd2q.m regarding the format.
[m n] = size(imag);
qr = imag(1,2:n);
qz = imag(2:m,1);
Int = imag(2:m,2:n);
delta_qz = qz(2)-qz(1);
delta_qr = qr(2)-qr(1);

%Set up the matrices whose elements are the values of q and phi for interpolation. 
q = q_range(1):q_delta:q_range(2);
phi = deg2rad(phi_range(1)):deg2rad(phi_delta):deg2rad(phi_range(2));
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

% Sum intensity along phi's for each q and create a column vector, I, as a
% function of q
I = sum(I,1)/length(phi);

% Create the output consisting of two column vectors.
result = [q' I'];

end