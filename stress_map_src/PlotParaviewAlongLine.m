clc
clear all
close all



fig = figure;
fig.Units  = 'centimeters';
fig.Position(3) = 9;
fig.Position(4) = 9;
set(fig.Children,'FontName','Times','FontSize',10);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space
fig.PaperPositionMode   = 'auto';
axis square
%set(gcf, 'color', 'none'); set(gca, 'color', 'none');

% 
% Refi= '../SelfStressInf2/FEMBounLineRef.csv' ;
% Refi = importdata(Refi)
% Refi = [Refi.data(:,4) Refi.data(:,20) Refi.data(:,22) Refi.data(:,21) Refi.data(:,17)];
% MakeMaps1(Refi(:,2)/1e6, Refi(:,3)/1e6, Refi(:,4)/1e6, Refi(:,5)/1e6, Refi(:,1)/10000, 1)

% Refi= '../SelfStressInf2/FEMBounLineUnif.csv' ;
% Refi = importdata(Refi)
% Refi = [Refi.data(:,4) Refi.data(:,20) Refi.data(:,22) Refi.data(:,21) Refi.data(:,17)];
% MakeMaps1(Refi(:,2)/1e6, Refi(:,3)/1e6, Refi(:,4)/1e6, Refi(:,5)/1e6, Refi(:,1)/10000, 1)


Refi= '../SelfStressInf2/FEMBoundLineRefWeygand.csv' ;
Refi = importdata(Refi)
Refi = [Refi.data(:,4) Refi.data(:,20) Refi.data(:,22) Refi.data(:,21) Refi.data(:,17)];
MakeMaps1(Refi(:,2)/1e6, Refi(:,3)/1e6, Refi(:,4)/1e6, Refi(:,5)/1e6, Refi(:,1)/10000, 1)



function MakeMaps1(SigmaXX,SigmaYY,SigmaXY,SigmaZZ, y,G)
 
    if G==1
    plot(y,SigmaXX,'k-','linewidth', 2)
    hold on
    plot(y,SigmaYY,'b-','linewidth', 2)
    plot(y,SigmaXY,'r-','linewidth', 2)
    plot(y,SigmaZZ,'g-','linewidth', 2)
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}')
    end
    if G==2
    hold on
    plot(y,SigmaXX,'.k')
    plot(y,SigmaYY,'b.')
    plot(y,SigmaXY,'r.')
    plot(y,SigmaZZ,'g.')
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}','FEM')
    end
    
    %hold off
    legend box off
    grid on
    set(gca,'FontWeight','bold');
    set(gca,'linewidth',1)
    xlabel('X (\mum)'); ylabel('\sigma (MPa)'); 

    xlim([-1 1])
    %print('-dpdf','-r600', 'AlongLineBoundary')
    print('-djpeg','-r600', 'FEMLineRefWeygand')  
end