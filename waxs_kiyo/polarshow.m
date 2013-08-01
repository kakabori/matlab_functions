function [q,phi,result] = polarshow(imag, q_range, phi_range, colorscale, q_delta, phi_delta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Get Global variables
global X_cen Y_cen; global X_Lambda Spec_to_Phos;

if (nargin<5)
    q_delta = 0.0025;
end

if (nargin<6)
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

result = interp2(imag,Y,X,'spline');

A = find(X > 1024);
B = find(X < 1);
C = find(Y > 1024);
D = find(Y < 1);
result(A) = 0;
result(B) = 0;
result(C) = 0;
result(D) = 0;

phi = rad2deg(phi);

colormap gray; imagesc(q,phi,result,colorscale);
axis xy;
set(gca,'xtick',min(q):(max(q)-min(q))/10:max(q));
set(gca,'xticklabel',min(q):(max(q)-min(q))/10:max(q));
set(gca,'xminortick','on');
set(gca,'ytick',min(phi):10:max(phi));
set(gca,'yticklabel',min(phi):10:max(phi));
set(gca,'yminortick','on');
set(gca,'tickdir','out');

end

