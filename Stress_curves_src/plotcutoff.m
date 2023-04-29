clc
close all
clear all

C = [0 5 10 25 50 75 100 200 300 400 500 600 1200]
Y = [194.78 188.421 188.267 187.672 186.379 185.203 184.037...
    178.159 178.501 176.788 180.206 182.011 180.648]

plotforpaper
hold on


plot(C,Y,'ko-', 'linewidth', 2)
xlabel('Cutoff distance (  )')% ($\AA$)', 'interpreter', 'latex')
ylabel('Yield stress (MPa)')
grid on
        
function plotforpaper
    fig = figure;
    fig.Units  = 'centimeters';
    fig.Position(3) = 8;
    fig.Position(4) = 8;
    set(fig.Children,'FontName','Times','FontSize',12);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space?
    fig.PaperPositionMode   = 'auto';
    %axis square?
    %axis image
    %set(gcf, 'color', 'none'); set(gca, 'color', 'none');
    set(gca,'FontWeight','bold');
    xlim([-20 1300])
    %ylim([140 220])
end