%function color_log2(range, nlines)
%
%Example1: color_log2([0 2000]);
%Example2:	color_log2([0 2000], 20);
%
%This function displays 2D image in gray log scale in a "filled contour map" way (contour "region")
%without drawing any actual contour lines. The program imitates contourf function. 
%range sets the color range, min and max. nlines determines 
%the number of levels of filled contour regions. The default value is 10.

function color_log2(range, nlines)

if (nargin < 2)
    nlines = 10;
end

%x sets the total number of contour regions. 
%Decrease stepsize to have more levels of contour colors.
stepsize = 1 / nlines;
x = [0:stepsize:1];
%g sets the gray color in logarithmic scale.
g = log(1:256)/log(256);
index = 256.^(x);
rindex = round(index);
N=length(rindex);
%This for loop fills in each range of color levels with the same grayness to imitate contourf function.
for k = 1:(N-1)
    g((rindex(k)+1):rindex(k+1)) = g(rindex(k+1));
end
cmap = [g' g' g'];
colormap(cmap)
caxis(range);
end
