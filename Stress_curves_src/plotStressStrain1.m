clc
close all
clear all

[file,path] = uigetfile({'*.csv';'*.txt'});
if isequal(file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(path,file)]);
end

 E = importdata(fullfile(path,file));
% 
% % %  "Time" "Real Time" 
% % %  "Sigma[11]($3)" "Sigma[22]($4)" "Sigma[33]($5)" "Sigma[12]($6)" "Sigma[13]($7)" "Sigma[23]($8)" 
% % % %  "E[11]($9)" "E[22]($10)" "E[33]($11)" "E[12]($12)" "E[13]($13)" "E[23]($14)" 
% % %  "Ep[11]($15)" "Ep[22]($16)" "Ep[33]($17)" "Ep[12]($18)" "Ep[13]($19)" "Ep[23]($20)"

%%%%%% plot strain stress and plastic strain in same graph SIGEPS.csv
if strcmp(file,'SIGEPS.csv')
    plotforpaper
    yyaxis left
    plot(E.data(:,11), E.data(:,5),'-.','linewidth', 2)
    xlabel('Strain')
    ylabel('Stress (MPa)')
    yyaxis right
    plot(E.data(:,11), E.data(:,17),'linewidth', 2)
    ylabel('Plastic Strain')
    grid on
end


% % % %   step[1] time[2]  Sigma33NC[3]   Sigma33C[4]  Str33NC[5]  Str33C[6]
% % % %   Pstr[7]  ElmerStrain[8] ElmerStress[9]

%%%%%% plot strain stress and plastic strain in same graph FEMLOAD.txt
if strcmp(file,'FEMLOAD.txt')
    plotforpaper
    yyaxis left
    plot(E.data(:,5),E.data(:,4),'-', 'linewidth', 2)
    xlabel('Strain')
    ylabel('Stress (MPa)')
    yyaxis right
    plot(E.data(:,5),E.data(:,7),'-', 'linewidth', 2)
    ylabel('Plastic Strain')
    grid on
end




function plotforpaper
    fig = figure;
    fig.Units  = 'centimeters';
    fig.Position(3) = 9;
    fig.Position(4) = 9;
    set(fig.Children,'FontName','Times','FontSize',12);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space?
    fig.PaperPositionMode   = 'auto';
    %axis square?
    axis image
    %set(gcf, 'color', 'none'); set(gca, 'color', 'none');
    set(gca,'FontWeight','bold');
end

%set(gcf, 'color', 'none'); set(gca, 'color', 'none');

