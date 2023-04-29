clc
close all

clear all

%A = importdata('DATA/FEMStress037.csv');
A = importdata('DATA/1micromDistance/FEMStressuniform.csv');

%%%% Figure settings
fig = figure;
fig.Units  = 'centimeters';
fig.Position(3) = 9;
fig.Position(4) = 9;
set(fig.Children,'FontName','Times','FontSize',10);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space
fig.PaperPositionMode   = 'auto';
axis square
%set(gcf, 'color', 'none'); set(gca, 'color', 'none');
plot(A.data(:,3)/10000,A.data(:,6)*1e-6,'r-','LineWidth',1) % XX
hold on
plot(A.data(:,3)/10000,A.data(:,9)*1e-6,'b-','LineWidth',1) % YY
plot(A.data(:,3)/10000,A.data(:,7)*1e-6,'k-','LineWidth',1) % XY
plot(A.data(:,3)/10000,A.data(:,11)*1e-6,'g-','LineWidth',1)% ZZ 
plot(A.data(:,3)/10000,A.data(:,8)*1e-6,'c-','LineWidth',1) % XZ
plot(A.data(:,3)/10000,A.data(:,10)*1e-6,'m-','LineWidth',1)% YZ

% plot(A.data(:,4)/10000,A.data(:,17)*1e-6,'r-','LineWidth',1) % XX
% hold on
% plot(A.data(:,4)/10000,A.data(:,20)*1e-6,'b-','LineWidth',1) % YY
% plot(A.data(:,4)/10000,A.data(:,18)*1e-6,'k-','LineWidth',1) % XY
% plot(A.data(:,4)/10000,A.data(:,22)*1e-6,'g-','LineWidth',1)% ZZ 
% plot(A.data(:,4)/10000,A.data(:,19)*1e-6,'c-','LineWidth',1) % XZ
% plot(A.data(:,4)/10000,A.data(:,21)*1e-6,'m-','LineWidth',1)% YZ
xlabel('{\it z} (µm)'); %,'fontweight','bold','fontsize',10)
ylabel('Stress (MPa)')
grid on

%B = importdata('DATA/GW3.dat')
B = importdata('DATA/1micromDistance/FivelData.txt')
plot(B.data(:,3),B.data(:,5),'k.','MarkerSize',9) % XY
plot(B.data(:,3),B.data(:,4),'r.','MarkerSize',9) % XX
plot(B.data(:,3),B.data(:,7),'b.','MarkerSize',9) % YY
plot(B.data(:,3),B.data(:,9),'g.','MarkerSize',9) % ZZ
plot(B.data(:,3),B.data(:,6),'c.','MarkerSize',9) % XZ
plot(B.data(:,3),B.data(:,8),'m.','MarkerSize',9) % YZ


legend('GW \sigma_{xx}','GW \sigma_{yy}','GW \sigma_{xy}','GW \sigma_{zz}','GW \sigma_{xz}','GW \sigma_{yz}','Sim', 'location', 'southeast')
legend box off
%grid on
xlim([0 2])
set(gca,'FontWeight','bold');
set(gca,'linewidth',1)

%print('-dpng','-r1000', 'FivelPlot037')
%print('-dpdf','-r600', 'Stresstest')
%export_fig('Trial.eps','-transparent');
