function a = slurp_rotate( image, angle )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a = slurp(image,'c');
a = imrotate(a, angle);

end

