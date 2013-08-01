%function result = sector_ccd(imag, q_range, phi_range, q_delta, phi_delta)
%   sect_water=sector_ccd(water,[0.6 2.1],[44.5 45.5],0.0025,0.1);
%
%This function will return intensity as a function of q along a constant value
%of phi on a CCD image. The output will contain two column vectors, 
%q and I. 
%The program will first transform the input 2D intensity map in 
%Cartensian coordinates into a map in polar coordinates for a region 
%specified by q_range and phi_range in steps of q_delta and phi_delta. 
%Then, the intensity will be averaged over a range specified by phi_range. 
%Data points will be created at every q_delta for a range specified 
%by q_range. 
%If no averaging is desired, phi_range(1) should be set equal to phi_range(2).
%
%Caveat: Since the input image will be interpolated, a user should pay
%attention to the actual resolution of the x-ray image and avoid
%interpolating many points between two real data points. The resolution in
%phi is sensitive to the position on the CCD detector; very close
%to the beam, phi resolution is very coarse. The resolution in q is almost
%constant over a whole CCD image for a typical sample-to-detector distance.
%   
%   Inputs  
%       imag -> N*M x-ray scattering image in Cartesian coordinates.

%       q_range -> create plot for q_range(1) < q < q_range(2).
%       Units are inverse Angstroms.
%       
%       phi_range -> will integrate intesity for phi_range(1) < phi < phi_range(2).
%       Angles are in degrees and measured counter-clockwise 
%       from the positive equator (right side of beam).
%                    
%       q_delta -> Resolution in q for the plot. Units are inverse Angstroms.
%       Should not be much
%       smaller than the resolution of the input image itself. 
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
%   Global variables -> X_cen, Y_cen , X_Lambda, Spec_to_Phos
%
%   KA 1/30/2012. 



function result = sector_ccd(imag, q_range, phi_range, q_delta, phi_delta)

% Get Global variables
global X_cen Y_cen; global X_Lambda Spec_to_Phos;

% Set resolution of interpolation in q
if (nargin<4)
    q_delta = 0.1;
end

% Set resolution of interpolation in phi
if (nargin<5)   
    phi_delta = 0.1;
end

q = q_range(1):q_delta:q_range(2);
phi = deg2rad(phi_range(1)):deg2rad(phi_delta):deg2rad(phi_range(2));
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

X = X_cen - y;
Y = Y_cen + x;

I = interp2(imag,Y,X,'spline');

A = find(X > 1024);
B = find(X < 1);
C = find(Y > 1024);
D = find(Y < 1);
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