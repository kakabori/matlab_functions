function [ output_args ] = saveqr(x,y,filename)
%saveqr This function will save the output from qrplot (2D plot)
%as a text file so that a user can import the data into Origin for graphing
%purposes.
%   x should be a row vector that will be used as the x coordinate. y is a
%   row vector that will be used as the y coordinate. This should be the
%   output from qrplot. filename must be a string, i.e. enclosed
%   by ' and '.
y_temp=y';
x_temp=x';
A=[x_temp y_temp];
save(filename,'A','-ascii');
end

