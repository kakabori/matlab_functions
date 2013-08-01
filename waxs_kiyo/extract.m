function [qr, qz, Int] = extract(imag)
%Extract qr and qz labels, and the intensity matrix according to the format
%specified by KA.
[m n] = size(imag);
qr = imag(1,2:n);
qz = imag(2:m,1);
Int = imag(2:m,2:n);

end

