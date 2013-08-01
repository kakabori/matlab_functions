function I_0 = apply_absorp_correct(imag, t, mu, alpha_d)
%This program applys the absorption correction to the input image. It assumes
%that the image is a fixed-angle data. The path length is calculated exactly.
%
%The program calculates the absorption of x-rays dues to the sample geometry
%for a given sample thickness, attenuation length, and angle of incidence. 
%Then, it calculates the absorption correction, which is simply the inverse of
%the absorption, multiplies the data by the correction, and returns the 
%corrected image.
%
%I_0 = apply_absorp_correct(imag, t, mu, alpha_d);
%correctedImage = abs_corr(uncorrectedImage, 10, 2300, 0.5);
%
%In the above example, the sample is 10 micron thick water layer at 0.5 degree
%incident angle, with 10.5 keV x-ray.
%Refer to the following website for a x-ray attenuation length for a given
%material: http://henke.lbl.gov/optical_constants/atten2.html
%
%INPUT
% imag -> CCD experimental data
% t -> thickness of the sample (microns)
% mu -> attenuation length of X-ray in the sample (microns)
% alpha_d -> angle of incidence (degrees)
%
%REQUIRED GLOBAL VARIABLES
% X_cen, Y_cen, Spec_to_Phos
%
%OUTPUT
% I_0 -> matrix representing the absorption-corrected image



global X_cen Y_cen; global Spec_to_Phos;

a = deg2rad(alpha_d);

[m n] = size(imag);

px = 1:n;
x = px - Y_cen;
len_x = length(x);

pz = 1:m;  
z = X_cen - pz;
len_z = length(z);

x2 = repmat(x,[len_z,1]);
z2 = repmat(z',[1,len_x]);

phi = atan(z2./x2);
th = 0.5 * atan(sqrt( (x2.^2+z2.^2)/Spec_to_Phos^2 ));

A = find(-sin(a)*cos(2*th) + cos(a)*sin(2*th).*sin(phi) < 0);
%unphysical result, that is, a negative path length of a photon after a 
%scattering event, will be stored in the variable A.
g = 1/sin(a) + 1./( -sin(a)*cos(2*th) + cos(a)*sin(2*th).*sin(phi) );

I_meas = (mu/t) * (1 - exp(-g*t/mu)) ./ g;
I_meas(A) = 0;%elements with unphysical result will be assigned a value zero.

I_0 = imag ./ I_meas;
I_0(A) = 0;
I_0(X_cen,Y_cen)=0;

end

