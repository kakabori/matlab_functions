function [ output_args ] = qshow(struct,color_scale)
%Last updated: 7/16/2012 by KA
%qshow(struct,color_scale)
%
%Use    qshow(q_water,[0 3000]);
%
%This function plots 2D intensity map in
%q-space. Use cshow if plotting a CCD image is desired (i.e. you want
%labels to be pixel numbers just like in tview program). It is particularly
%useful when one wants to plot an image that is transformed into q-space
%from CCD space because the dimension of a matrix is not 1024 by 1024
%anymore. The function that does the transformation (transform_ccd2q.m) outputs
%a structure containing qr and qz vectors as well as the transformed image. For each Matlab
%function, consult help in Matlab.

imagesc(struct.qr,struct.qz,struct.Int,color_scale);
colormap(gray); axis image; axis xy; %axis([1.3 2.0 0 0.5]);

%Notation. [a:b:c] means tick marks will be placed between a and c with an
%increment of b. For example, [-1:0.5:1] will generate marks at -1, -0.5,
%0, 0.5, and 1.
set(gca,'xtick',min(struct.qr):0.1:max(struct.qr));
set(gca,'ytick',min(struct.qz):0.1:max(struct.qz));
set(gca,'xticklabel',min(struct.qr):0.1:max(struct.qr));
set(gca,'yticklabel',min(struct.qz):0.1:max(struct.qz));
set(gca,'tickdir','out');
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
set(gca,'FontSize',10);
xlabel(' q_r (Å^{-1})','FontSize',12);
ylabel(' q_z (Å^{-1})','FontSize',12);
end

