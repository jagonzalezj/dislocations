clc
close all
clear all

%         Azul gordo             rojo mandarina          verde olivo
C = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880], ...
    [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};
%      Amarillo naranja            violeta              Azul cielo claro       rojo hierro sangre

plotforpaper
hold on

for i = 1:4

    [file,path] = uigetfile({'*.dat'});
    E = importdata(fullfile(path,file));


    if i == 1 | i == 2
        file2 = 'SIGEPS.csv'
    else
        file2 = 'FEMLOAD.txt'
    end

    D = importdata(fullfile(path,file2));
    lower = min(length(D.data(:,2)),length( E.data(:,2)))

    
    if strcmp(file2,'FEMLOAD.txt')
        
        %   step[1] time[2]  Sigma33NC[3]   Sigma33C[4]  Str33NC[5]  Str33C[6]
        %   Pstr[7]  ElmerStrain[8] ElmerStress[9]
%         yyaxis right
%         plot(D.data(:,2),D.data(:,7),'--','color',C{i}, 'linewidth', 1);
        
        %yyaxis left
        plot(D.data(:,5),D.data(:,4),'--','color',C{i}, 'linewidth', 1);


    end
    
    
    if strcmp(file2,'SIGEPS.csv')
        %   1         2           3           4           5           6
        % "Time" "Real Time" "Sigma[11]" "Sigma[22]" "Sigma[33]" "Sigma[12]" 
        %       7         8          9       10      11       12      13
        % "Sigma[13]" "Sigma[23]" "E[11]" "E[22]" "E[33](" "E[12](" "E[13]" 
        %   14        15      16      17        18       19        20    21
        % "E[23]" "Ep[11]" "Ep[22]" "Ep[33]" "Ep[12]" "Ep[13]" "Ep[23]" "sequence"
        
%         yyaxis right
%         plot(D.data(:,2), D.data(:,17),'--','color',C{i}, 'linewidth', 1)
        
        %yyaxis left
        plot(D.data(:,11), D.data(:,5),'--','color',C{i}, 'linewidth', 1)

        
        
    end
    
%     ylabel('Plastic Strain')
%     
%     % 1.#STEP 2.#EL_ENERGY (J) 3.#CORE_ENERGY (J) 4.#SF_ENERGY (J) 
%     % 5.#T_ENERGY (J) 6.#SUR_AREA (A^2) 7.#DIS_LENGTH (A)
%     yyaxis left
%     plot(D.data(1:lower,2), E.data(1:lower,7),'-','color',C{i},'linewidth', 1)
%     ylabel('Disloc length (Ang)')

end

%xlabel('Time (ns)')


    
legend('NumC0','NumC600','ElNumC0','ElNumC600')




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