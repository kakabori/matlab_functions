%Last updated: 8/10/2011 by KA
function [ output_args ] = cshow(image,color_scale)
%cshow This function opens an empty figure and plots a 1024-by-1024
%intensity map in CCD-space. Use qshow if plotting an image in q-space is
%desired. Use show if you want labels to be the actual matrix row and
%column number. cshow is particularly useful when one wants to plot an
%image that is equivalent to images displayed in tview program. show is
%perhaps more useful for data analysis. For each Matlab function, consult
%help in Matlab.

px=[1:1024];%px sets the label for x axis. 
py=[1:1024];
figure;
colormap gray;
imagesc(px,py,image,color_scale);
axis image;
%axis xy;

%Notation. [a:b:c] means tick marks will be placed between a and c with an
%increment of b. For example, [-1:0.5:1] will generate marks at -1, -0.5,
%0, 0.5, and 1.
set(gca,'xtick',[200:200:1000]);%decide where tick marks appear for x axis
set(gca,'ytick',[25:200:825]);%decide where tick marks appear for y axis
set(gca,'xticklabel',[200:200:1000]);%label for x axis.
set(gca,'yticklabel',[1000:-200:200]);%label for y axis.
set(gca,'tickdir','out');%tick mark direction: out or in
set(gca,'xminortick','on');
set(gca,'yminortick','on');
xlabel(' p_x ','FontSize',16);
ylabel(' p_z ','FontSize',16);
end
