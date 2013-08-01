function AddQrQzLabels(qr, qz)
%addQrQzLabels(qr, qz) 
%
%This function adds qr and qz labels to the current figure.
%A user must supply qr and qz row vectors which contain qr and qz values
%at which tick marks and labels will be placed. These vectors must be
%monotonically increasing order, but the increments do not have to be 
%the same throughout the vectors.
%
%addQrQzLabels(0:0.1:1, [0.2 0.3 0.8]);
%will create labels at qr = 0, 0.1, 0.2, ..., 1 and qz = 0.2, 0.3, and 0.8.

global lambda beamX beamZ sDist

qz = fliplr(qz);
theta = asin(lambda * qr / 4 /pi);
px = beamX + sDist * tan(2*theta);
theta = asin(lambda * qz / 4 /pi);
pz = beamZ - sDist * tan(2*theta);
set(gca, 'FontName', 'Arial', 'FontSize', 12);
set(gca, 'xtick', px);
set(gca, 'ytick', pz);
set(gca, 'xticklabel', qr);
set(gca, 'yticklabel', qz);
set(gca, 'TickDir', 'out');
%xlabel('q_r (Å^{-1})');



end

