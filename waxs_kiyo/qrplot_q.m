%Last updated: 2/16/2012
function result = qrplot_q(imag, qz_range)
%qrplot_ka 
%   

[m n] = size(imag);
qr = imag(1,2:n);
qz = imag(2:m,1);
Int = imag(2:m,2:n);
delta_qr = qr(2)-qr(1);

A = find(qz >= qz_range(1) & qz < (qz_range(1)+delta_qr));
B = find(qz >= qz_range(2) & qz < (qz_range(2)+delta_qr));

Int2 = Int(A:B,:);
Int3 = mean(Int2,1);

result = [qr,Int3];
fprintf('Intensity was integrated from qz=%g A^-1 to qz=%g A^-1,\n',qz(A),qz(B));
fprintf('centered at %g\n',(qz(A)+qz(B))/2);
plot(qr',Int3');
end

