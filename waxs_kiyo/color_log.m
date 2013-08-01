%function color_log(range)
%
%Use    color_log([0 11000]);
%
%This function displays 2D image in gray log scale.

function color_log(range)

g = log(1:256)/log(256);
cmap = [g' g' g'];
colormap(cmap);
caxis(range);
end
