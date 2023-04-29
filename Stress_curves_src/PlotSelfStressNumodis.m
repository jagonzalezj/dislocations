%%%------------------------------------------------------------------------
% Plot Self Stress and FEM stress in a plane guiven it coordiates points
%  -How to Use:
%     1. Run a ELmerNumodis timestep and with Paraview Export a planar
%        spreadsheet with all parameters (Pos Displ Stress and Strain)
%        with file name FEMStressSurface.csv
%     2. Runing this script will save the node position in 1d array(pos.txt) 
%        as demanded for numodis inquire the self stress. The  option 
%        SaveData = 1 otherwise pos.txt will not be saved and only the FEM
%        stress will be plotted.
%     3. Run the same ELmerNumodis algorithm, like that the couplig will
%        read the pos.txt file and export the self stress in a file named
%        SelfStress.txt.
%     4. Running this script one more time with option Comparative = 1 will
%        read the SelfStress.txt file and plot the values as we all compare
%        them (FEM-SELF)
%
%   Two plotting options are available (one with 4 components for image 
%   stress comparison at perpendicular planar faces as Hirt and Lothe image
%   calculation quasi2d) and another with the six components for plotting
%   at the free surfaces.

%%%------------------------------------------------------------------------

clc
close all
clear all
format long g

dist = 1;   % distance to surface: if [0.2 microm] => [dist = 2] 

% discretization in fivel benchmark
%disc1 = 50;  % fem discretization in x
%disc2 = 50;  % fem discretization in y

% discretization in inf edge benchmark
disc1 = 150;  % fem discretization in x
disc2 = 150;  % fem discretization in y

%%%%%-----------------------   DATA COMING FROM PARAVIEW    ---------------

% load the paraview spreadsheet data .csv
%A = importdata(['FEMStressSurface0_',num2str(dist),'microm.csv']);
%A = importdata('../Fivel2/FEMStressBottom.csv');

A = importdata('FEMStressBottom.csv');

% slice the postion data
pos = [A.data(:,1) A.data(:,2) A.data(:,3) A.data(:,15:20)];

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
[Sxx, Syy, Szz, Sxy, Sxz, Syz] = SplitStressComp(pos(:,4:9)*1e-6, 2, disc1, disc2);

% Separate data position  and components units conversions (A -> microm)
Xm = reshape(pos(:,c(1))*1e-4, [disc2, disc1]);
Ym = reshape(pos(:,c(2))*1e-4, [disc2, disc1]);

% Send to plot FEM Stress
  %MakeMaps2(Sxx, Syy, Szz, Sxy, Sxz, Syz, Xm, Ym, 1, [-10 10]);
MakeMaps3(cat(6,Sxx, Syy, Szz, Sxy, Sxz, Syz), Xm, Ym, 1, [-20 20]);

%Txt save the data
 if isfile(['pos0_',num2str(dist),'microm.txt']) == 0; %SaveData == 1
     disp '*** Saving point coordinate data for numodis Couplig...';
     disp '*** Saving mode: total#  x1 y1 z1 x2 y2 z2 x3 y3 z3 ...';
     point4numodis = reshape(pos(:,1:3)', [1 3*length(pos(:,1))]);
     fid1=fopen('pos.txt','w');
     fprintf(fid1, '%d \n', length(point4numodis));
     fprintf(fid1, '%f ', point4numodis');
     fclose(fid1);
 end

%%%%%-----------------------   DATA COMING FROM NUMODIS    ---------------

%%%Compare to Self Stress
direccion = 'SelfStress.txt' ;
%direccion1= '../SelfStressInf2/WeygandStressUnif.txt' ;
if isfile(direccion) == 1 % & Comparative == 1
    disp 'Pos file exit runing comparative script';
    
    Self = load(direccion);
    Self = AdaptStress(Self);
    [Slfxx, Slfyy, Slfzz, Slfxy, Slfxz, Slfyz] = SplitStressComp(Self, 1, disc1, disc2);
    
%     Wey = load(direccion1);
%     Wey = AdaptStress(Wey);
%     [Weyxx, Weyyy, Weyzz, Weyxy, Weyxz, Weyyz] = SplitStressComp(Wey, 1, disc1, disc2);
    
    %plot Self Stress
  %  MakeMaps2(Slfxx, Slfyy, Slfzz, Slfxy, Slfxz, Slfyz, Xm, Ym, 2, [-10 10]);
    %MakeMaps3(cat(6,Slfxx, Slfyy, Slfzz, Slfxy, Slfxz, Slfyz), Xm, Ym, 2, [-20 20])
   MakeMaps3(cat(6,Slfxx, Slfyy, Slfzz, Slfxy, Slfxz, Slfyz), Xm, Ym, 2, [-20 20])
    
    %plot Weygand Stress
    %MakeMaps2(Weyxx, Weyyy, Weyzz, Weyxy, Weyxz, Weyyz, Xm, Ym, 3, [-20 20]);
    %MakeMaps3(cat(6,Weyxx, Weyyy, Weyzz, Weyxy, Weyxz, Weyyz), Xm, Ym, 3, [-20 20])
    %MakeMaps3(cat(6,Weyxx, Weyyy, Weyzz, Weyxy, Weyxz, Weyyz), Xm, Ym, 4, [-20 20])

    %Plot FEMstress + SelfStress + Weygand
    %MakeMaps2(Sxx+Slfxx, Syy+Slfyy, Szz+Slfzz, Sxy+Slfxy, Sxz+Slfxz, Syz+Slfyz, Xm, Ym, 4, [-10 10]);
    %MakeMaps2(Sxx-Slfxx-Weyxx, Syy-Slfyy-Weyyy, Szz-Slfzz-Weyzz, Sxy-Slfxy-Weyxy, Sxz-Slfxz-Weyxz, Syz-Slfyz+Weyyz, Xm, Ym, 4, [-20 20]);
    %MakeMaps2(Slfxx+Weyxx, Slfyy+Weyyy, Slfzz+Weyzz, Slfxy+Weyxy, Slfxz+Weyxz, Slfyz+Weyyz, Xm, Ym, 4, [-20 20]);
    %MakeMaps2(Sxx+Weyxx, Syy+Weyyy, Szz+Weyzz, Sxy+Weyxy, Sxz+Weyxz, Syz+Weyyz, Xm, Ym, 4, [-20 20]);
    MakeMaps3(cat(6,Sxx+Slfxx, Syy+Slfyy, Szz+Slfzz, Sxy+Slfxy, Sxz+Slfxz, Syz+Slfyz), Xm, Ym, 4, [-20 20]);

    %MakeMaps3(cat(6,Sxx+Slfxx+Weyxx, Syy+Slfyy+Weyyy, Szz+Slfzz+Weyzz, Sxy+Slfxy+Weyxy, Sxz+Slfxz+Weyxz, Syz+Slfyz+Weyyz), Xm, Ym, 4, [-20 20]);
end

disp("Finished")

clearvars -except fictX fictY


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


%%%-----------------------------------------------------------------------
function MakeMaps2(Sxx, Syy, Szz, Sxy, Sxz, Syz, Xm, Ym, G, param3)

    function plotmap1 % this function specify the caxis  :)
        if param3(1) == 0

        else
            caxis([param3(1) param3(2)])
        end
    end


    h1=figure(G);

    subplot(2,3,1);
    %imagesc(Syy)
    %imagesc(Xm(1,:), Ym(:,1), Sxx); set(gca,'YDir','normal')

    contourf(Xm,-Ym, Sxx, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Sxx (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1
    xlim([-1 1])
    ylim([-1 1])

    subplot(2,3,2); 
    contourf(Xm,Ym, Syy, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Syy (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1
        xlim([-1 1])
    ylim([-1 1])

    subplot(2,3,3); 
    contourf(Xm,Ym, Szz, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Szz (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1
        xlim([-1 1])
    ylim([-1 1])

    subplot(2,3,4); 
    contourf(Xm,Ym, Sxy, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Sxy (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1
        xlim([-1 1])
    ylim([-1 1])

    subplot(2,3,5); 
    contourf(Xm,Ym, Sxz, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Sxz (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1
        xlim([-1 1])
    ylim([-1 1])

    subplot(2,3,6); 
    contourf(Xm,Ym, Syz, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Syz (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1
        xlim([-1 1])
    ylim([-1 1])

    colormap(jet)
    set(h1,'position',[10 10 1000 1000]);

    names = {'FEM Stress', 'Self Stress', 'Weygand Stress','FEM + Self Stress', 'SelfStress + ImageStress',...
        'Airy + SelfStress + ImageStress'};

    suptitle(names(G));
    
    %print('-dpng','-r1000', G)
end
%%%-----------------------------------------------------------------------

%%%-----------------------------------------------------------------------
function MakeMaps3 (S, Xm, Ym, G, param3)
nombres = ['Sxx'; 'Syy' ;'Szz'; 'Sxy'; 'Sxz'; 'Syz'];
maping = {'FEM_Stress', 'SelfNumEdge', 'Weygand_Stress','FEM + Self_Stress', 'SelfStress + ImageStress',...
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
        %pcolor(S(:,:,o))
        %shading interp

        %colormap jet
        %caxis([-20 20])
       
        contourf(Xm,Ym, S(:,:,o), 150, 'ShowText', 'off', 'LineColor', 'none')

        %title(['',num2str(maping{G}),'',num2str(nombres(o,:)),'']); 
        %colorbar
        %xlabel('X (µm)'); ylabel('Y (µm)'); 
        axis image
        axis off

        %set(gca,'FontWeight','bold');
        plotmap1
        colormap(jet)
        
        xlim([-0.5 0.5])
        ylim([-0.5 0.5])
       
        print('-dpng','-r1000', ['',num2str(maping{G}),'-',num2str(nombres(o,:)),''])
    end
end
%%%-----------------------------------------------------------------------