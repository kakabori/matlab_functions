function [ output_args ] = qcontour(struct,color_scale)
%Last updated: 7/16/2012 by KA
%qcontour(struct,color_scale)
%
%Use    qcontour(q_water,[0 3000]);
%
%This function plots 2D color filled contour map in q-space. The function that does the 
%transformation (transform_ccd2q.m) outputs a structure containing qr and qz vectors as well as 
%the transformed image. For each Matlab function, consult help in Matlab.

[cout, H] = contourf(struct.qr,struct.qz,struct.Int,color_scale);
colormap(gray(256)); axis image; %axis xy; %axis([1.3 2.0 0 0.5]);
colorbarf(cout, H);

%Notation. [a:b:c] means tick marks will be placed between a and c with an
%increment of b. For example, [-1:0.5:1] will generate marks at -1, -0.5,
%0, 0.5, and 1.
set(gca,'xtick',min(struct.qr):0.1:max(struct.qr));
set(gca,'ytick',min(struct.qz):0.1:max(struct.qz));
set(gca,'yticklabel',min(struct.qz):0.1:max(struct.qz));
set(gca,'xticklabel',min(struct.qr):0.1:max(struct.qr));
set(gca,'tickdir','out');
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
set(gca,'FontSize',10);
xlabel(' q_r (Å^{-1})','FontSize',12);
ylabel(' q_z (Å^{-1})','FontSize',12);
end
