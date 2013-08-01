%Last updated: 2/16/2012
function result = qzplot_q(imag, qr_range)
%qzplot_ka 
%   

[m n] = size(imag);
qr = imag(1,2:n);
qz = imag(2:m,1);
Int = imag(2:m,2:n);
delta_qz = qz(2)-qz(1);

A = find(qr >= qr_range(1) & qr < (qr_range(1)+delta_qz));
B = find(qr >= qr_range(2) & qr < (qr_range(2)+delta_qz));

Int2 = Int(:,A:B);
Int3 = mean(Int2,2);

result = [qz Int3];
fprintf('Intensity was integrated from qr=%g A^-1 to qr=%g A^-1,\n',qr(A),qr(B));
fprintf('centered at %g\n',(qz(A)+qz(B))/2);
plot(qz,Int3);
end
