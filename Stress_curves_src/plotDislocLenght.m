clc
close all
clear all

C = {'b', 'r', 'g','k'};
plotforpaper
hold on

for i = 1:4

    [file,path] = uigetfile({'*.dat'});

    E = importdata(fullfile(path,file));

    % 1.#STEP 2.#EL_ENERGY (J) 3.#CORE_ENERGY (J) 4.#SF_ENERGY (J) 
    % 5.#T_ENERGY (J) 6.#SUR_AREA (A^2) 7.#DIS_LENGTH (A)

    plot(E.data(:,1), E.data(:,7),'-','linewidth', 2)

end


xlabel('Steps')
ylabel('Disloc length (Ang)')
    
legend('NumW0','NumW600','ElNumW0','ElNumW600')




function plotforpaper
    fig = figure;
    fig.Units  = 'centimeters';
    fig.Position(3) = 9;
    fig.Position(4) = 9;
    set(fig.Children,'FontName','Times','FontSize',12);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space?
    fig.PaperPositionMode   = 'auto';
    %axis square?
    %axis image
    %set(gcf, 'color', 'none'); set(gca, 'color', 'none');
    set(gca,'FontWeight','bold');
    %xlim([1.3e-3 2.5e-3])
    %ylim([140 220])
end