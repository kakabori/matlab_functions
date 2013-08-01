function [ output_args ] = saveqz(x,y,filename)
%saveqr This function will save the output from qzplot (2D plot) as a text
%file so that a user can import the data into Origin for graphing purposes.
%   x should be a row vector that will be used as the x coordinate. y is a
%   column vector that will be used as the y coordinate. This should be the
%   output from qzplot. filename must be a string, i.e. enclosed by two 's.
x_temp=x';
A=[x_temp y];
save(filename,'A','-ascii');
end

