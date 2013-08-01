function contourplot(image)
%% function contourplot
%
%

%global cout H vector;

flip = flipud(image);
%vector = [-1000 0 50*exp(0:0.6:6)];
vector = [-1000 0 50 91 170 300 550 1000 1800 3300 6100 11000];
figure;
[cout,H,cf]=contourf(flip,vector);
axis image;
%vector2=[0 20 54 150 400 1100 3000 8000];
%vector2=[-100 0 50*exp(0:0.6:6)];
vector2 = [-100 0 50 91 170 300 550 1000 1800 3300 6100 11000];
color_nonlinear(vector2);
%colormap(gray(256));
set(H,'edgecolor','none');
axis ([200,800,100,900]);
%colorbarf(cout,H);
qlabel_contourf;

