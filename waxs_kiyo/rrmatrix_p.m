function rrmatrix_p( filein, p_range, eta_range, d_eta, ccdtype )
%% RingRockMatrix(rrmatrix) creates omega-eta matrix for mosaicity analysis.
% rrmatrix_p(filein, p_range, eta_range, d_eta, ccdtype)
%
% rrmatrix_p('cmu_DOPC_Tat_10to1_2012.txt', [446 450], [-10 10], 0.1, 'r');
% 
% The three outputs are out_matrix.txt, out_somega.txt, and out_seta.txt.
% omega is the angle of incidence and eta is the angle defined from the meridian.
%
% The program will overwrite the existing files with the above file names
% in the current working directory, so be sure to change the names of the 
% files if you want to keep them.
%
% The required global variables are X_cen and Y_cen.
%
%
% filein -> The name of the plain text file that contains the information
%           about the files that will be analyzed. The first column must be the
%           file names. The second column must be the angle of incidence. They must
%           be seprated by TAB.
%
% p_range -> [pmin pmax], the range over which intensity on a ring will be averaged.
%            The range is measured from the beam, i.e., (pmin+pmax)/2 is the radius in 
%            pixels measured from the beam, along which a ring will be produced.
% 
% eta_range -> [min max], the range over which the ring curves will be produced
% 
% d_eta -> step size for eta. 0.1 works well.
%
% ccdtype -> 'r' for data taken at CMU with Rigaku CCD. 'c' for data taken 
%            at CHESS. When 'r' is supplied, the program will use slurpNagle function.
%            When 'c' is supplied, the program will use slurp function.

%% Get file names and angles
fid = fopen(filein);
A = textscan(fid, '%s %f');
fclose(fid);
filename = A{1};
omega = A{2}';
omega2 = [NaN omega];

%% ccdtype (CMU Rigaku or Chess)
if ccdtype == 'r'
    image = slurpNagle(filename{1});
    int_image = ringplot(image,p_range,eta_range,d_eta);
    intensity = zeros(length(int_image),(length(filename)));
    intensity(:,1) = int_image(:,2);

    for n = 2:length(filename)
        image = slurpNagle(filename{n});
        int_image = ringplot(image,p_range,eta_range,d_eta);
        intensity(:,n) = int_image(:,2);
    end

elseif ccdtype == 'c'
    image = slurp(filename{1},'c');
    image = imrotate(image, -90); %option for chess04a dmpc5 data
    int_image = ringplot(image,p_range,eta_range,d_eta);
    intensity = zeros(length(int_image),(length(filename)));
    intensity(:,1) = int_image(:,2);

    for n = 2:length(filename)
        image = slurp(filename{n},'c');
        image = imrotate(image, -90); %option for chess04a dmpc5 data
        int_image = ringplot(image,p_range,eta_range,d_eta);
        intensity(:,n) = int_image(:,2);
    end

else
    fprintf('\nChoose the correct ccd type (Rigaku or Chess):\n');
    fprintf('''r'' for CMU data       ''c'' for CHESS data\n');
end

%% Save matrix.txt
eta = int_image(:,1);
ring = [eta intensity];
ring = [omega2; ring];

global X_cen Y_cen;

fid = fopen('out_matrix.txt', 'w');
fprintf(fid, 'The second line, beginning with NaN, is the values of omega. It should go to the comment line in OriginPro.\n');
fprintf(fid, 'This data was produced with rrmatrix_p(''%s'',[%d %d],[%.1f %.1f],%.1f,''%s'');\n',filein,p_range(1),p_range(2),eta_range(1),eta_range(2),d_eta,ccdtype);
fprintf(fid, 'X_cen = %d; Y_cen = %d;\n', X_cen, Y_cen);
fprintf(fid, '      eta,      intensity\n');
fclose(fid);
dlmwrite('out_matrix.txt', ring, 'delimiter', ',', 'precision', '%9.2f', '-append');

%% Save somega.txt
somega = zeros(length(eta)*length(omega), 3);

for k = 1:length(omega)
    somega(((k-1)*length(eta)+1):(k*length(eta)),1) = eta;
    somega(((k-1)*length(eta)+1):(k*length(eta)),2) = omega(k);
    somega(((k-1)*length(eta)+1):(k*length(eta)),3) = intensity(:,k);
end

fid = fopen('out_somega.txt', 'w');
fprintf(fid, 'This data was produced with rrmatrix_p(''%s'',[%d %d],[%.1f %.1f],%.1f,''%s'');\n',filein,p_range(1),p_range(2),eta_range(1),eta_range(2),d_eta,ccdtype);
fprintf(fid, 'X_cen = %d; Y_cen = %d;\n', X_cen, Y_cen);
fprintf(fid, '%9s,%9s,  %s\n','eta', 'omega', 'intensity');
fclose(fid);
dlmwrite('out_somega.txt', somega, 'delimiter', ',', '-append', 'precision', '%9.2f');
    
%% Save seta.txt
seta = zeros(length(eta)*length(omega), 3);

for k = 1:length(eta)
    seta(((k-1)*length(omega)+1):(k*length(omega)),1) = eta(k);
    seta(((k-1)*length(omega)+1):(k*length(omega)),2) = omega;
    seta(((k-1)*length(omega)+1):(k*length(omega)),3) = intensity(k,:);
end

fid = fopen('out_seta.txt', 'w');
fprintf(fid, 'This data was produced with rrmatrix_p(''%s'',[%d %d],[%.1f %.1f],%.1f,''%s'');\n',filein,p_range(1),p_range(2),eta_range(1),eta_range(2),d_eta,ccdtype);
fprintf(fid, 'X_cen = %d; Y_cen = %d;\n', X_cen, Y_cen);
fprintf(fid, '%9s,%9s,  %s\n','eta','omega','intensity');
fclose(fid);
dlmwrite('out_seta.txt', seta, 'delimiter', ',', '-append', 'precision', '%9.2f');