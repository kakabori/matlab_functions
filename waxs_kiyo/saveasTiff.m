function [ output_args ] = saveasTiff( q_image , filename )
%% saveasTiff saveasTiff(q_water);
%   q_image -> N*M x-ray scattering image in q-space Cartesian coordinates.
%           The input image must be an output of transform_ccd2q function so
%           that it follows the specific format defined by KA. See
%           transform_ccd2q.m for the details of the format.

%% Extract labels
%Extract qr and qz labels, and the intensity matrix according to the format
%specified by KA.
[m n] = size(q_image);
qr = q_image(1,2:n);
qz = q_image(2:m,1);
Int = q_image(2:m,2:n);

%%
Int = uint16(Int);
imwrite(Int,filename);

end

