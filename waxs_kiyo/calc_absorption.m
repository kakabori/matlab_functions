% What is this function for?

function [ I_meas ] = absorption(t, mu, alpha_d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global X_cen Y_cen; global Spec_to_Phos;

a=deg2rad(alpha_d);

px = 1:1024;
x = px - Y_cen;
len_x = length(x);

pz = 1:1024;  
z = X_cen - pz;
len_z = length(z);

x2 = repmat(x,[len_z,1]);
z2 = repmat(z',[1,len_x]);

phi = atan(z2./x2);
th = 0.5 * atan(sqrt( (x2.^2+z2.^2)/Spec_to_Phos^2 ));

A = find(-sin(a)*cos(2*th) + cos(a)*sin(2*th).*sin(phi) < 0);
g = 1/sin(a) + 1./( -sin(a)*cos(2*th) + cos(a)*sin(2*th).*sin(phi) );

I_meas = (mu/t) * (1 - exp(-g*t/mu)) ./ g;

I_meas = 1 ./ I_meas;
I_meas(A) = 0;
I_meas(X_cen,Y_cen)=0;


%figure; colormap gray;
%imagesc(px,pz,I_meas,colorscale);

%imagesc(px,px,Int);
%axis image;axis xy;


end

