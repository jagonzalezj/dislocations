%%%%% this one uses the system of references rotated 
%                  ^
%representation    | = Z;  --> = Y;  .=Z
clc 
close all
clear all

discX = 60;
discY = 25;
l = 0.1;  % microm

A = importdata('../Weygand_study2NonPBC/ResultSlab/NoWeygand60.csv');

%%%%%  <-------------POSTION------------->|<----STRESES--->
pos = [A.data(:,1) A.data(:,2) A.data(:,3) A.data(:,21:26)];

% 1. Sort by horizontal lines in Z [depend of orientation]
pos = sortrows(pos, [3]);


% 2.1  IF ONLY ONE REFINEMENT
%    Apply first filter [ the consecutives discY points are replaced every 
%    by their comon mean value enssuring a right organization of points
for i = 1:discX;
    pos((discY*i)-(discY-1):discY*i,3) = mean(pos((discY*i)-(discY-1):discY*i,3));
end

% 2.2 IF REFINEMENT IN THE TWO DIMENTIONS
%    Apply first filter [ the consecutives discY points are replaced every 
%    by their comon mean value enssuring a right organization of points
% for i = 1:discY;
%     pos((discX*i)-(discX-1):discX*i,2) = mean(pos((discX*i)-(discX-1):discX*i,2));
% end

% 3. Now sort firt by 3 and then 2 column (data fully organized)
pos = sortrows(pos, [3,2]); 

% 4. Reshape of columns to matrixes to contournf them
%    A---> micrometer     Pa----> MPa
  Y = reshape(pos(:,3),[discY,discX])'*1e-4;
  X = reshape(pos(:,2),[discY,discX])'*1e-4;
Sxx = reshape(pos(:,4),[discY,discX])'*1e-6;
Sxy = reshape(pos(:,5),[discY,discX])';
Sxz = reshape(pos(:,6),[discY,discX])';
Syy = reshape(pos(:,7),[discY,discX])'*1e-6;
Syz = reshape(pos(:,8),[discY,discX])'*1e-6;
Szz = reshape(pos(:,9),[discY,discX])'*1e-6;

%MakeMaps1(Syz, Syy, Szz, Sxx, X, Y, 1, [-0 0],l);
MakeMaps1(Syy, Szz, Syz, Sxx, X-0.5, Y, 1, [0 0],l);

%MakeMaps2(Sxx, Syy, Szz, Sxy, Sxz, Syz, X, Y, 1, [-0 0]);
%MakeMaps3(cat(6,Sxx, Syy, Szz, Sxy, Sxz, Syz), X, Y, 1, [-0 0]);





%%%-----------------------------------------------------------------------

function MakeMaps1(SigmaXX,SigmaYY,SigmaXY,SigmaZZ, x, y, G, param3, l)

    function plotmap1 % this function specify the caxis  :)
        if param3(1) == 0
        else
            caxis([param3(1) param3(2)])
        end
    end

    h1=figure(G);
    subplot(2,2,1); 
    contourf(x,y, SigmaZZ, 150, 'ShowText', 'off', 'LineColor', 'none')
    %imagesc(x(1,:), y(:,1), SigmaZZ); set(gca,'YDir','normal')
    %xline(0,'--','linewidth',2)
    %hold on
    %plot(-l,0,'k.','markersize', 10)
    title('\sigma_{zz} (MPa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');colorbar;
    plotmap1

    subplot(2,2,2); 
    contourf(x,y, SigmaXX, 150, 'ShowText', 'off', 'LineColor', 'none')
    %imagesc(x(1,:), y(:,1),SigmaXX); set(gca,'YDir','normal')
    %xline(0,'--','linewidth',2)
    %hold on
    %plot(-l,0,'k.','markersize', 10)
    title('\sigma_{xx} (MPa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');colorbar;
    plotmap1

    subplot(2,2,3);
    contourf(x,y, SigmaYY, 150, 'ShowText', 'off', 'LineColor', 'none')
    %imagesc(x(1,:), y(:,1), SigmaYY); set(gca,'YDir','normal')
    %xline(0,'--','linewidth',2)
    %hold on
    %plot(-l,0,'k.','markersize', 10)
    title('\sigma_{yy} (MPa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');
    colorbar
    plotmap1

    subplot(2,2,4); 
    contourf(x,y, SigmaXY, 150, 'ShowText', 'off', 'LineColor', 'none')
    %imagesc(x(1,:), y(:,1), SigmaXY); set(gca,'YDir','normal')
    %xline(0,'--','linewidth',2)
    %hold on
    %plot(-l,0,'k.','markersize', 10)
    title('\sigma_{xy} (MPa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold'); colorbar;
    plotmap1

    set(h1,'position',[10 10 1000 1000]); colormap(jet);
    
    names = {'SelfStress', 'ImageStress','Airy function', 'SelfStress + ImageStress',...
            'Airy + SelfStress + ImageStress', 'Airy+ ImageStress'};
%     suptitle(names(G));
    
%     if G == 6    
%     print('-dpng','-r600', 'AiryImage')
%     end
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

    subplot(2,3,2); 
    contourf(Xm,Ym, Syy, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Syy (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1

    subplot(2,3,3); 
    contourf(Xm,Ym, Szz, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Szz (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1

    subplot(2,3,4); 
    contourf(Xm,Ym, Sxy, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Sxy (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1

    subplot(2,3,5); 
    contourf(Xm,Ym, Sxz, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Sxz (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1

    subplot(2,3,6); 
    contourf(Xm,Ym, Syz, 150, 'ShowText', 'off', 'LineColor', 'none')
    title('Syz (Mpa)'); colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis image
    set(gca,'FontWeight','bold');
    plotmap1

    colormap(jet)
    set(h1,'position',[10 10 1000 1000]);

    names = {'FEM Stress', 'Self Stress', 'Weygand Stress','FEM + Self Stress+ Weygang', 'SelfStress + ImageStress',...
        'Airy + SelfStress + ImageStress'};

    suptitle(names(G));
    
    %print('-dpng','-r600', G)
end
%%%-----------------------------------------------------------------------

%%%-----------------------------------------------------------------------
function MakeMaps3 (S, Xm, Ym, G, param3)
nombres = ['Sxx'; 'Syy' ;'Szz'; 'Sxy'; 'Sxz'; 'Syz'];
maping = {'FEM_Stress', 'SelfNumEdge', 'Weygand_Stress','FEM + Self_Stress+ Weygang', 'SelfStress + ImageStress',...
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
       
        %print('-dpng','-r1000', ['Devinc_',num2str(maping{G}),'-',num2str(nombres(o,:)),''])
    end
end
%%%-----------------------------------------------------------------------
