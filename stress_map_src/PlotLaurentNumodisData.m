clc
close all
clear all


% discretization in inf edge benchmark
disc1 = 128;  % fem discretization in x
disc2 = 128;  % fem discretization in y


A = importdata('../laurent/screw/H&L.csv');


% slice the postion data
pos = [A.data(:,15) A.data(:,16) A.data(:,17) A.data(:,8:13)];

% detect the not Normal plane [1 3] means 2 is the normal
c = [0 0];
l= 0;
for j = 1:3
    if length(unique(pos(:,j))) ~= 1;
       l = l+1;
       c(l) = j; 
    end 
end

% create ficticious vectors (finish this) vector index; number of elements
fictX = min(round(pos(:,c(1)))):(max(round(pos(:,c(1))))-...
            min(round(pos(:,c(1)))))/(disc1-1):max(round(pos(:,c(1))));
        
fictY = min(round(pos(:,c(2)))):(max(round(pos(:,c(2))))-...
            min(round(pos(:,c(2)))))/(disc2-1):max(round(pos(:,c(2))));

% update row data postion with fict data (finish this) what columns to
% take
pos(:,c(1)) = raw2fict(pos(:,c(1)), fictX);
pos(:,c(2)) = raw2fict(pos(:,c(2)), fictY);

%  Sort Position matrix by row consecutivelly (Finish this) uses the same
%  coluns of previoys sentence
pos  = sortrows(pos, [c(1),c(2)]);

% Separate data components and convert FEM stress to MPa ojooooooo aquiiii
[Sxx, Syy, Szz, Sxy, Sxz, Syz] = SplitStressComp(pos(:,4:9), 1, disc1, disc2);

% Separate data position  and components units conversions (A -> microm)
Xm = reshape(pos(:,c(1))*1e-4, [disc2, disc1]);
Ym = reshape(pos(:,c(2))*1e-4, [disc2, disc1]);

% Send to plot FEM Stress
%  MakeMaps2(Sxx, Syy, Szz, Sxy, Sxz, Syz, Xm, Ym, 1, [-20 20]);
MakeMaps3(cat(6,Sxx, Syy, Szz, Sxy, Sxz, Syz), Xm, Ym, 1, [-20 20]);

%Txt save the data
%  if isfile(['pos0_',num2str(dist),'microm.txt']) == 0; %SaveData == 1
%      disp '*** Saving point coordinate data for numodis Couplig...';
%      disp '*** Saving mode: total#  x1 y1 z1 x2 y2 z2 x3 y3 z3 ...';
%      point4numodis = reshape(pos(:,1:3)', [1 3*length(pos(:,1))]);
%      fid1=fopen('pos.txt','w');
%      fprintf(fid1, '%d \n', length(point4numodis));
%      fprintf(fid1, '%f ', point4numodis');
%      fclose(fid1);
%  end


%%%%%-----------------------   DATA COMING FROM NUMODIS    ---------------

%%%Compare to Self Stress
%direccion = '../StrightEdgeScrewSlefStress/SelfStress.txt' ;
%direccion1= '../SelfStressInf2/WeygandStressUnif.txt' ;
%if isfile(direccion) == 1 % & Comparative == 1
%    disp 'Pos file exit runing comparative script';
    
%    Self = load(direccion);
%    Self = AdaptStress(Self);
%    [Slfxx, Slfyy, Slfzz, Slfxy, Slfxz, Slfyz] = SplitStressComp(Self, 1, disc1, disc2);
    
%     Wey = load(direccion1);
%     Wey = AdaptStress(Wey);
%     [Weyxx, Weyyy, Weyzz, Weyxy, Weyxz, Weyyz] = SplitStressComp(Wey, 1, disc1, disc2);
   

    %plot Self Stress
    %MakeMaps2(Slfxx, Slfyy, Slfzz, Slfxy, Slfxz, Slfyz, Xm, Ym, 2, [-20 20]);
    %MakeMaps3(cat(6,Slfxx, Slfyy, Slfzz, Slfxy, Slfxz, Slfyz), Xm, Ym, 2, [-20 20])
%    MakeMaps3(cat(6,Slfxx, Slfyy, Slfzz, Slfxy, Slfxz, Slfyz), Xm, Ym, 2, [-20 20])
    
    %plot Weygand Stress
    %MakeMaps2(Weyxx, Weyyy, Weyzz, Weyxy, Weyxz, Weyyz, Xm, Ym, 3, [-20 20]);
    %MakeMaps3(cat(6,Weyxx, Weyyy, Weyzz, Weyxy, Weyxz, Weyyz), Xm, Ym, 3, [-20 20])
    %MakeMaps3(cat(6,Weyxx, Weyyy, Weyzz, Weyxy, Weyxz, Weyyz), Xm, Ym, 4, [-20 20])

    %Plot FEMstress + SelfStress + Weygand
    %MakeMaps2(Sxx+Slfxx, Syy+Slfyy, Szz+Slfzz, Sxy+Slfxy, Sxz+Slfxz, Syz+Slfyz, Xm, Ym, 5, [0 0]);
    %MakeMaps2(Sxx-Slfxx-Weyxx, Syy-Slfyy-Weyyy, Szz-Slfzz-Weyzz, Sxy-Slfxy-Weyxy, Sxz-Slfxz-Weyxz, Syz-Slfyz+Weyyz, Xm, Ym, 4, [-20 20]);
    %MakeMaps2(Slfxx+Weyxx, Slfyy+Weyyy, Slfzz+Weyzz, Slfxy+Weyxy, Slfxz+Weyxz, Slfyz+Weyyz, Xm, Ym, 4, [-20 20]);
    %MakeMaps2(Sxx+Weyxx, Syy+Weyyy, Szz+Weyzz, Sxy+Weyxy, Sxz+Weyxz, Syz+Weyyz, Xm, Ym, 4, [-20 20]);

    %MakeMaps3(cat(6,Sxx+Slfxx+Weyxx, Syy+Slfyy+Weyyy, Szz+Slfzz+Weyzz, Sxy+Slfxy+Weyxy, Sxz+Slfxz+Weyxz, Syz+Slfyz+Weyyz), Xm, Ym, 4, [-20 20]);
%end

disp("Finished")

%clearvars -except fictX fictY




%%%-----------------------------------------------------------------------
function MakeMaps3 (S, Xm, Ym, G, param3)
nombres = ['Sxx'; 'Syy' ;'Szz'; 'Sxy'; 'Sxz'; 'Syz'];
maping = {'FEM_Stress', 'NumEdge', 'Weygand_Stress','FEM + Self_Stress+ Weygang', 'SelfStress + ImageStress',...
        'Airy + SelfStress + ImageStress'};

    function plotmap1 % this function specify the caxis  :)
        if param3(1) == 0

            else
                caxis([param3(1) param3(2)])
        end
    end

    for o=1:6
        nombres(o,:);
        fig = figure;
        fig.Units  = 'centimeters';
        fig.Position(3) = 6.0;
        fig.Position(4) = 6.0;
        set(fig.Children,'FontName','Times','FontSize',10);
        %set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space
        set(gca, 'Position', get(gca, 'OuterPosition')); %Yes this remove  the white space
        fig.PaperPositionMode   = 'auto';
               
        %imagesc(Slfyy)
        pcolor(S(:,:,o))
        shading interp
        %colormap jet
        %caxis([-20 20])      
        
        %contourf(Xm,Ym, S(:,:,o), 150, 'ShowText', 'off', 'LineColor', 'none')
        %title(['',num2str(maping{G}),'',num2str(nombres(o,:)),'']); 
        %colorbar
        %xlabel('X (µm)'); ylabel('Y (µm)'); 
        axis image
        axis off

        %set(gca,'FontWeight','bold');
        plotmap1
        colormap(jet)
       
        print('-dpng','-r1000', ['H&L_',num2str(maping{G}),'-',num2str(nombres(o,:)),''])
    end
end
%%%-----------------------------------------------------------------------




%%%-----------------------------------------------------------------------
% Adapt the real raw pos data from paraview to a one artificially created
%%%-----------------------------------------------------------------------
function rawvect = raw2fict(rawvect, fictvect)
    for k = 1: length(rawvect)
        [~,idx]=min(abs(fictvect-rawvect(k)));
        rawvect(k) = fictvect(idx);
    end
end

%%%-----------------------------------------------------------------------
% reshape the six component stress tensor (in column format) to a
% meshgrided format to be plotted. 2 ways can be used Numodis format and
% ElmerSpreadsheet format
%%%-----------------------------------------------------------------------
function [S1,S2,S3,S4,S5,S6] = SplitStressComp(tobesplit, modo, disc1, disc2)
    switch modo
        case 1 % numodis format
            S1 = reshape(tobesplit(:,1), [disc2,disc1]);
            S2 = reshape(tobesplit(:,2), [disc2,disc1]);
            S3 = reshape(tobesplit(:,3), [disc2,disc1]);
            S4 = reshape(tobesplit(:,4), [disc2,disc1]);
            S5 = reshape(tobesplit(:,5), [disc2,disc1]);
            S6 = reshape(tobesplit(:,6), [disc2,disc1]);
        case 2
            S1 = reshape(tobesplit(:,1), [disc2,disc1]);
            S2 = reshape(tobesplit(:,4), [disc2,disc1]);
            S3 = reshape(tobesplit(:,6), [disc2,disc1]);
            S4 = reshape(tobesplit(:,2), [disc2,disc1]);
            S5 = reshape(tobesplit(:,3), [disc2,disc1]);
            S6 = reshape(tobesplit(:,5), [disc2,disc1]);
    end
end

%%%-----------------------------------------------------------------------
% Adapt the 1d arrar of the stress components coming from Numodsi coupling
% to a six component column array
%%%-----------------------------------------------------------------------
function S = AdaptStress(Stress)
% Separate the six components of the SelfStress S(n,6)
m = 0;
S = [];
    for i = 1:length(Stress)/6
        m= m + 1;
        S(i,1) = Stress(6*m-5);
        S(i,2) = Stress(6*m-4);
        S(i,3) = Stress(6*m-3);
        S(i,4) = Stress(6*m-2);
        S(i,5) = Stress(6*m-1);
        S(i,6) = Stress(6*m);
    end
end
%%%-----------------------------------------------------------------------




