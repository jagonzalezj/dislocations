clc
close all
clear all


Coord = importdata('pos.txt'); % import one file data (numodis format)
Coord = Coord(2:end);         % avoid first number (total points)

% create the 3d vectors to plot.
for i = 1: length(Coord)/3
    X(i)=[Coord(3*i - 2)]/10000;
    Y(i)=[Coord(3*i - 1)]/10000;
    Z(i)=[Coord(3*i)]/10000;
end

% Optional : plot the plane where stress is requested
plot3(X,Y,Z,"k."); xlabel("X(1)"); ylabel("Y(2)"); zlabel("Z(3)"); grid on; axis square



%Here perape the pltting functions
repX = unique(X);
repY = unique(Y);
repZ = unique(Z);

Stress = load('SelfStress.txt');

if length(repX)== 1
    disp('ploting normal to YZ')
    [Xm, Ym] = meshgrid(repY, repZ);
    al = repY;
    be = repZ;
    
elseif length(repY) == 1
    disp('ploting normal to XZ')
    [Xm, Ym] = meshgrid(repX, repZ);    
    al = repX;
    be = repZ;
    
elseif length(repZ)==1
    disp('ploting normal to XY')
    [Xm, Ym] = meshgrid(repX, repY);
    al = repX; 
    be = repY;
end
    
m = 0;

    for i = 1 : length(al)
        for j = 1 : length(be)
            m= m + 1;
            Sxx(i,j) = Stress(6*m-5); 
            Syy(i,j) = Stress(6*m-4);
            Szz(i,j) = Stress(6*m-3);
            Sxy(i,j) = Stress(6*m-2);
            Sxz(i,j) = Stress(6*m-1);
            Syz(i,j) = Stress(6*m);
        end
    end

    
    %%%% ploting 
    
    h1=figure(2);
    subplot(3,2,1); contourf(Xm,Ym, Sxx, 20, 'ShowText', 'on')
    title('Sxx (Mpa): Inf. edge'); %colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis square
    
    subplot(3,2,2); contourf(Xm,Ym, Syy, 30, 'ShowText', 'on')
    title('Syy (Mpa): Inf. edge'); %colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis square
    
    subplot(3,2,3); contourf(Xm,Ym, Szz, 30, 'ShowText', 'on')
    title('Szz (Mpa): Inf. edge'); %colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis square
    
    subplot(3,2,4); contourf(Xm,Ym, Sxy, 30, 'ShowText', 'on')
    title('Sxy (Mpa): Inf. edge'); %colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis square
    
    subplot(3,2,5); contourf(Xm,Ym, Sxz, 30, 'ShowText', 'on')
    title('Sxz (Mpa): Inf. edge'); %colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis square
       
    subplot(3,2,6); contourf(Xm,Ym, Syz, 30, 'ShowText', 'on')
    title('Syz (Mpa): Inf. edge'); %colorbar
    xlabel('X (\mum)'); ylabel('Y (\mum)'); axis square
    colormap(jet)
    set(h1,'position',[10 10 1000 1000]);
    
    