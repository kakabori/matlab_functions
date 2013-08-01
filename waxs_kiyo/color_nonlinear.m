function color_nonlinear( x )

Nx = length(x);
clim = [min(x) max(x)];
dx = min(diff(x));
y = clim(1):dx:clim(2);
for k=1:Nx-1, y(y>x(k) & y<=x(k+1)) = x(k+1); end
cmap = colormap(gray(Nx));
cmap2 = [...
interp1(x(:),cmap(:,1),y(:)) ...
interp1(x(:),cmap(:,2),y(:)) ...
interp1(x(:),cmap(:,3),y(:)) ...
];
colormap(cmap2)
caxis(clim)
