function qlabel_contourf
%% function qlabel
%
%

%% Parse global variables
global X_cen Y_cen X_Lambda Spec_to_Phos

%Assuming 1024 by 1024 image,
xcen = 1024 - X_cen + 1;

%% Set up labels
qr=[-2.0:0.1:2.0];
qz2=[0:0.1:2];
qr2=[-2.0:0.2:2.0];

%% 
%ypixels = -Spec_to_Phos * tan( 2* asin( X_Lambda*qz/(4*pi) ) ) + X_cen;
ypixels = Spec_to_Phos * tan( 2* asin( X_Lambda*qz2/(4*pi) ) ) + xcen;
xpixels = Spec_to_Phos * tan( 2* asin( X_Lambda*qr/(4*pi) ) ) + Y_cen;
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
set(gca,'xtick',xpixels);
set(gca,'ytick',ypixels);
set(gca,'FontSize',10);
set(gca,'yticklabel',qz2);
set(gca,'xticklabel',abs(qr));
set(gca,'TickDir','out');
h1=xlabel('$\mathrm{q_r(\AA^{-1})}$','FontSize',14);
h2=ylabel('$\mathrm{q_z(\AA^{-1})}$','FontSize',14);
set(h1,'Interpreter','LaTex');
set(h2,'Interpreter','LaTex');
pos1=get(h1,'pos');
pos2=get(h2,'pos');
set(h1,'pos',pos1+[0 0 0]);
set(h2,'pos',pos2+[-10 0 0]);


%ylabel(' q_z (Å^{-1})','FontSize',12);
