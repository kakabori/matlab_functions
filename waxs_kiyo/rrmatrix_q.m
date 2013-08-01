function rrmatrix_q( filein, q_range, phi_range, delta_phi, ccdtype )
%% RingRockMatrix(rrmatrix) creates omega-eta matrix for mosaicity analysis.
%
% The three outputs are out_matrix.txt, out_somega.txt, and out_seta.txt.
%
% The program will overwrite the existing files with the above file names
% in the current working directory, so be sure to change the names of the 
% files if you want to keep them.
%
% filein -> The name of the plain text file that contains the information
% about the files that will be analyzed. The first column must be the
% file names. The second column must be the angle of incidence. They must
% be seprated by TAB.
%
% q_range -> An input for integrate_annulus_ccd
% 
% phi_range -> An input for integrate_annulus_ccd function
% 
% delta_phi -> An input for integrate_annulus_ccd function
%

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
    int_image = integrate_annulus_ccd(image,q_range,phi_range,delta_phi,1,0.0025);
    intensity = zeros(length(int_image),(length(filename)));
    intensity(:,1) = int_image(:,2);

    for n = 2:length(filename)
        image = slurpNagle(filename{n});
        int_image = integrate_annulus_ccd(image,q_range,phi_range,delta_phi,1,0.0025);
        intensity(:,n) = int_image(:,2);
    end

elseif ccdtype == 'c'
    image = slurp(filename{1},'c');
    int_image = integrate_annulus_ccd(image,q_range,phi_range,delta_phi,1,0.0025);
    intensity = zeros(length(int_image),(length(filename)));
    intensity(:,1) = int_image(:,2);

    for n = 2:length(filename)
     image = slurp(filename{n},'c');
        int_image = integrate_annulus_ccd(image,q_range,phi_range,delta_phi,1,0.0025);
        intensity(:,n) = int_image(:,2);
    end

else
    fprintf('\nChoose the correct ccd type (Rigaku or Chess):\n');
    fprintf('''r'' for CMU data       ''c'' for Chess data\n');
end

%% Save matrix.txt
eta = 90 - int_image(:,1);
ring = [eta intensity];
ring = [omega2; ring];

fid = fopen('out_matrix.txt', 'w');
fprintf(fid, 'The second line, beginning with NaN, is the values of omega. It should go to the comment line in OriginPro.\n');
fprintf(fid, 'This data was produced with rrmatrix_q(''%s'',[%.3f %.3f],[%.1f %.1f],%.1f,''%s'');\n',filein,q_range(1),q_range(2),phi_range(1),phi_range(2),delta_phi,ccdtype);
fprintf(fid, '     eta,     intensity\n');
fclose(fid);
dlmwrite('out_matrix.txt', ring, 'delimiter', ',', 'precision', '%8.1f', '-append');

%% Save somega.txt
somega = zeros(length(eta)*length(omega), 3);

for k = 1:length(omega)
    somega(((k-1)*length(eta)+1):(k*length(eta)),1) = eta;
    somega(((k-1)*length(eta)+1):(k*length(eta)),2) = omega(k);
    somega(((k-1)*length(eta)+1):(k*length(eta)),3) = intensity(:,k);
end

fid = fopen('out_somega.txt', 'w');
fprintf(fid, '%8s,%8s,  %s\n','eta', 'omega', 'intensity');
fclose(fid);
dlmwrite('out_somega.txt', somega, 'delimiter', ',', '-append', 'precision', '%8.1f');
    
%% Save seta.txt
seta = zeros(length(eta)*length(omega), 3);

for k = 1:length(eta)
    seta(((k-1)*length(omega)+1):(k*length(omega)),1) = eta(k);
    seta(((k-1)*length(omega)+1):(k*length(omega)),2) = omega;
    seta(((k-1)*length(omega)+1):(k*length(omega)),3) = intensity(k,:);
end

fid = fopen('out_seta.txt', 'w');
fprintf(fid, '%8s,%8s,  %s\n','eta','omega','intensity');
fclose(fid);
dlmwrite('out_seta.txt', seta, 'delimiter', ',', '-append', 'precision', '%8.1f');

 
 
 