function [ output_args ] = image2d(imag,colorscale)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

colormap gray; imagesc(imag,colorscale);
axis xy;

end

