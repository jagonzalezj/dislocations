clc
close all

clear all

A = importdata('FEMStress1.csv');

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

plot(A.data(:,3)/10000,A.data(:,6)*1e-6,'m-','LineWidth',1)
hold on
plot(A.data(:,3)/10000,A.data(:,7)*1e-6,'g-','LineWidth',1)
plot(A.data(:,3)/10000,A.data(:,8)*1e-6,'c-','LineWidth',1)
plot(A.data(:,3)/10000,A.data(:,9)*1e-6,'r-','LineWidth',1)
plot(A.data(:,3)/10000,A.data(:,10)*1e-6,'k-','LineWidth',1)
plot(A.data(:,3)/10000,A.data(:,11)*1e-6,'b-','LineWidth',1)
xlabel('Z (µm)'); %,'fontweight','bold','fontsize',10)
ylabel('\sigma (MPa)')


B = importdata('marc.txt');
plot(B(:,3),B(:,9),'b.','LineWidth',2)
plot(B(:,3),B(:,4),'m.','LineWidth',2)
plot(B(:,3),B(:,5),'g.','LineWidth',2)
plot(B(:,3),B(:,6),'c.','LineWidth',2)
plot(B(:,3),B(:,7),'r.','LineWidth',2)
plot(B(:,3),B(:,8),'k.','LineWidth',2)


legend('\sigma_{xx}','\sigma_{xy}','\sigma_{xz}','\sigma_{yy}','\sigma_{yz}','\sigma_{zz}','GW method', 'location', 'southeast')
legend box off
%grid on
xlim([0 2])
set(gca,'FontWeight','bold');
set(gca,'linewidth',1)


%print('-dpdf','-r600', 'Stresstest')
%export_fig('Trial.eps','-transparent');
