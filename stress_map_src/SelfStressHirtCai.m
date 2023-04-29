clc
clear all
close all

% % Material properties: Copper (Numodis)
%E = 111216 ;            % Young Modulus  [Mpa]
psi = 0.324;             % Poisson ratio
b = 0.25526;              % Burgers vector [nm]
miu = 42000;             % [MPa]
E = 2*miu*(1+psi);
b = b * 0.001;     % nm    ==>  micrometer

l=0.1;



% Domain definition z perpendicular view direction Theory 
X = [-1:0.05:3];  % micrometers
Y = [-7.5:0.05:7.5];
[x,y]=meshgrid(X,Y);
%plot(x,y,'b.')
%%%%-----------------------------------------------------------------------
%%%%%%      Self Stress  Edge (Hirt and Lothe {Singular})          %%%%%%%%
%%%%-----------------------------------------------------------------------
% prim = (miu*b)/(2*pi*(1-psi));
% SigmaXX =  -prim * ( (y.*(3*x.^2 + y.^2))./ (x.^2 + y.^2).^2 );
% SigmaYY =   prim * ( (y.*(x.^2 - y.^2))./ (x.^2 + y.^2).^2 );
% SigmaXY =   prim * ( (x.*(x.^2 - y.^2))./ (x.^2 + y.^2).^2 );
% SigmaZZ =  -prim *(2*psi) * (y./(x.^2 + y.^2)); 
%%SigmaZZ = psi*(SigmaXX + SigmaYY);

%%%%%%     Self Stress  Screw (Hirt and Lothe {Singular})
% SigmaXZ = -(miu*b)/(2*pi) * (y./(x.^2 + y.^2));
% SigmaYZ = -(miu*b)/(2*pi) * (x./(x.^2 + y.^2));
%%%%-----------------------------------------------------------------------

% Domain definition z perpendicular view direction Theory 
X = [-1:0.05:0];  % micrometers
Y = [-1:0.05:1];
[x,y]=meshgrid(X,Y);
x=x+l;  % move the disloc to -l postion
%%%%-----------------------------------------------------------------------
%%%%%%           Self Stress  Edge (Cai {Non-Singular})            %%%%%%%%
%%%%-----------------------------------------------------------------------
a = 0.0003;
prim = (miu*b)/(2*pi*(1-psi));
rhoa = sqrt(a^2 + x.^2 + y.^2);
SigmaXX =  -prim * (y./rhoa.^2).*(1 + (2*(x.^2 + a^2)./rhoa.^2));
SigmaYY =   prim * (y./rhoa.^2).*(1 - (2*(y.^2 + a^2)./rhoa.^2));
SigmaXY =   prim * (x./rhoa.^2).*(1 - ((2*(y.^2))./rhoa.^2));
SigmaZZ =  -prim * (2*psi) .* (y./rhoa.^2).*(1 + (a^2./rhoa.^2));
%%%SigmaZZ = psi*(SigmaXX + SigmaYY)

%%%%%%     Self Stress  Screw (Cai {Non-Singular}
% SigmaXZ = -(miu*b)/(2*pi) * (y./rhoa.^2) .* (1 + (a^2/rhoa.^2));
% SigmaYZ = -(miu*b)/(2*pi) * (x./rhoa.^2) .* (1 + (a^2/rhoa.^2));
% SELF STRESS PLOT
x=x-l;  % acomodate the plot x axe
h1=figure(1);
subplot(2,2,1); 
contourf(x, y, SigmaZZ,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Self ZZ (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');
colorbar;
%caxis([-10 10])

subplot(2,2,2); 
contourf(x, y, SigmaXX,50, 'ShowText', 'off', 'LineColor', 'none')
title('Self XX (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');
colorbar; 
%caxis([-15 15])

subplot(2,2,3); 
contourf(x, y, SigmaYY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Self YY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold'); 
%caxis([-15 15])
colorbar

subplot(2,2,4); contourf(x, y, SigmaXY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Self XY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold'); 
colorbar;
%caxis([-15 15])

set(h1,'position',[10 10 1000 1000]); colormap(jet);
%print -dpng -r600 InfEdgHirtCu1micronSxx1


% Domain definition z perpendicular view direction Theory 
X = [0:0.05:1];  % micrometers
Y = [-1:0.05:1];
[x,y]=meshgrid(X,Y);
x=x-l;  % move the disloc to +l position
%%%%-----------------------------------------------------------------------
%%%%%%         IMAGE  Self Stress  Edge (Cai {Non-Singular})            %%%%%%%%
%%%%-----------------------------------------------------------------------
a = 0.0003;
prim = (miu*(-b))/(2*pi*(1-psi));
rhoa = sqrt(a^2 + x.^2 + y.^2);
ISigmaXX =  -prim * (y./rhoa.^2).*(1 + (2*(x.^2 + a^2)./rhoa.^2));
ISigmaYY =   prim * (y./rhoa.^2).*(1 - (2*(y.^2 + a^2)./rhoa.^2));
ISigmaXY =   prim * (x./rhoa.^2).*(1 - ((2*(y.^2))./rhoa.^2));
ISigmaZZ =  -prim * (2*psi) .* (y./rhoa.^2).*(1 + (a^2./rhoa.^2));
%%%SigmaZZ = psi*(SigmaXX + SigmaYY)

%%%%%%     Self Stress  Screw (Cai {Non-Singular}
% SigmaXZ = -(miu*b)/(2*pi) * (y./rhoa.^2) .* (1 + (a^2/rhoa.^2));
% SigmaYZ = -(miu*b)/(2*pi) * (x./rhoa.^2) .* (1 + (a^2/rhoa.^2));
% SELF STRESS PLOT
x=x+l;  % acomodate the plot x axe
h2=figure(2);
subplot(2,2,1); 
contourf(x, y, ISigmaZZ,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Im ZZ (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');colorbar;caxis([-10 10])

subplot(2,2,2); 
contourf(x, y, ISigmaXX,50, 'ShowText', 'off', 'LineColor', 'none')
title('Im XX (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');colorbar; caxis([-15 15])

subplot(2,2,3); 
contourf(x, y, ISigmaYY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Im YY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold'); caxis([-15 15])
colorbar

subplot(2,2,4); 
contourf(x, y, ISigmaXY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Im XY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold'); colorbar;caxis([-15 15])

set(h2,'position',[10 10 1000 1000]); colormap(jet);


x=x-max(X)/2;
h3=figure(3);
subplot(2,2,1); 
contourf(x, y, ISigmaZZ+SigmaZZ,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('\sigma_{ZZ} (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');colorbar;caxis([-10 10])
hold on; xline(0); hold off

subplot(2,2,2); 
contourf(x, y, ISigmaXX+SigmaXX,50, 'ShowText', 'off', 'LineColor', 'none')
title('\sigma_{XX} (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');colorbar; caxis([-15 15])
hold on; xline(0); hold off

subplot(2,2,3); 
contourf(x, y, ISigmaYY+SigmaYY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('\sigma_{YY} (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');colorbar; caxis([-15 15])
hold on; xline(0); hold off

subplot(2,2,4); 
contourf(x, y, ISigmaXY+SigmaXY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('\sigma_{XY} (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold'); colorbar;caxis([-15 15])
hold on; xline(0); hold off

set(h3,'position',[10 10 1000 1000]); colormap(jet);
%print -dpng -r1000 SurfaceFree

% 
% 
% 
% 
% 
% 
% 
% 
% 
%%%%----------------------------------------------------------------------
%%%%%%   Calculation of image stress (Theory from Hirt and Lothe)
%%%%%%   There is:
%%%%%%       Total_Stress 
%%%%%%       (Total_Stress) - Self_Stress = Image_Stress
%%%%-----------------------------------------------------------------------
%%% Acomodate the distance to surface to mathc with self stress plot
%%% the "0" is at the boundarie, so change it :)
X = [-1:0.05:0];  % micrometers
Y = [-1:0.05:1];
[x,y]=meshgrid(X,Y);
%x = x+1; % micrometers
%[x,y]=meshgrid(X,Y);
%%%%-----------------------------------------------------------------------
%%%%%%      Airy Edge (Hirt and Lothe {Singular}) 
%%%%-----------------------------------------------------------------------
l = 0.1;
r = sqrt ((x-l).^2 + y.^2);
I_SigmaXX = -prim * (4*l*x.*y)./(r.^6) .* (3.*(l-x).^2 - y.^2);
I_SigmaYY =  prim * (2*l)./(r.^6) .* ( 4*y.*(l-x).^3 + 6*x.*y.*(l-x).^2 + 4*(y.^3).*(l-x) - 2*x.*y.^3);
I_SigmaXY = -prim * (2*l)./(r.^6) .* ( (l-x).^4 + 2*x.*(l-x).^3 - 6*x.*(y.^2).*(l-x) - y.^4);
I_SigmaZZ =  prim * (8*l*psi)./(r.^6) .* (y.*(l-x).^3 + (l-x).*(y).^3); 
%%%%-----------------------------------------------------------------------


%--% PLOTIING Total_Stress 
%%%%-----------------------------------------------------------------------
h2=figure(4);
subplot(2,2,2); contourf(x, y, I_SigmaXX,20, 'ShowText', 'off', 'LineColor', 'none'); %hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Airy stress XX (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
set(gca,'FontWeight','bold');
axis equal; 
colorbar

subplot(2,2,3); contourf(x, y, I_SigmaYY,20, 'ShowText', 'off', 'LineColor', 'none'); %hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Airy stress YY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
set(gca,'FontWeight','bold');
axis equal; 
colorbar

subplot(2,2,1); contourf(x, y, I_SigmaZZ, 20, 'ShowText', 'off', 'LineColor', 'none'); %hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Airy stress ZZ (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
set(gca,'FontWeight','bold');
axis equal; 
colorbar

subplot(2,2,4); contourf(x, y, I_SigmaXY,20, 'ShowText', 'off', 'LineColor', 'none'); %hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Airy stress XY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mumm)'); 
set(gca,'FontWeight','bold');
axis equal; 
colorbar
colormap(jet)
set(h2,'position',[10 10 1000 1000]);

MakeMaps3(cat(4,I_SigmaZZ, I_SigmaXX, I_SigmaYY, I_SigmaXY), x, y, 1, [-20 20]);
%print -dpng -r1000 AiryEdg01Micron

% 
% 
% %--% PLOTIING(Total_Stress) - Self_Stress = Image_Stress
% %%%%-----------------------------------------------------------------------
% h3=figure(3);
% subplot(2,2,2); contourf(x, y, I_SigmaXX+SigmaXX,20, 'ShowText', 'off', 'LineColor', 'none');hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
% title('I_SigmaXX (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
% axis equal; 
% colorbar
% 
% subplot(2,2,3); contourf(x, y, I_SigmaYY+SigmaYY,20, 'ShowText', 'off', 'LineColor', 'none');hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
% title('I_SigmaYY (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
% axis equal; 
% colorbar
% 
% subplot(2,2,1); contourf(x, y, I_SigmaZZ+SigmaZZ, 20, 'ShowText', 'off', 'LineColor', 'none');hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
% title('I_SigmaZZ (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
% axis equal; 
% colorbar
% 
% % solve the names using latex structure
% subplot(2,2,4); contourf(x, y, I_SigmaXY+SigmaXY,20, 'ShowText', 'off', 'LineColor', 'none');hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
% title('I_SigmaXY (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mumm)'); 
% axis equal; 
% colorbar
% colormap(jet)
% set(h3,'position',[10 10 1000 1000]); 
% %%%%-----------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% %--% PLOTIING DIFFERENCES(Self_Stress-Num) - Self_Stress_Cai
% %%%%-----------------------------------------------------------------------
% % % h3=figure(4);
% % % subplot(2,2,1); contourf(x, y,SigmaZZ-Sxx, 20, 'ShowText', 'off');hold on; %plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
% % % title('I_SigmaZZ (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); axis square
% % % set(gca,'FontWeight','bold'); colorbar;
% % % caxis([-1 1])
% % % 
% % % subplot(2,2,2); contourf(x, y,SigmaXX-Syy,20, 'ShowText', 'off');hold on; %plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
% % % title('I_SigmaXX (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); axis square
% % % set(gca,'FontWeight','bold'); colorbar;
% % % caxis([-1 1])
% % % 
% % % subplot(2,2,3); contourf(x, y,SigmaYY-Szz,20, 'ShowText', 'off');hold on; %plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
% % % title('I_SigmaYY (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); axis square
% % % set(gca,'FontWeight','bold'); colorbar;
% % % caxis([-1 1])
% % % 
% % % subplot(2,2,4); contourf(x, y,SigmaXY-Syz,20, 'ShowText', 'off');hold on; %plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
% % % title('I_SigmaXY (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mumm)'); axis square
% % % set(gca,'FontWeight','bold'); colorbar;
% % % caxis([-1 1])
% % % colormap(jet)
% % % set(h3,'position',[10 10 1000 1000]);
% % % print -dpng -r600 InfEdgDiff4Micron
% 
% %%%%-----------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Stress of the image dislocation in the side of the real  
X = [-1:0.05:0];  % micrometers
Y = [-1:0.05:1];
[x,y]=meshgrid(X,Y);%
x=x-l;  % move the disloc to +l position
%%%%-----------------------------------------------------------------------
%%%%%%         IMAGE  Self Stress  Edge (Cai {Non-Singular})            %%%%%%%%
%%%%-----------------------------------------------------------------------
a = 0.0003;
prim = (miu*(-b))/(2*pi*(1-psi));
rhoa = sqrt(a^2 + x.^2 + y.^2);
ISigmaXX =  -prim * (y./rhoa.^2).*(1 + (2*(x.^2 + a^2)./rhoa.^2));
ISigmaYY =   prim * (y./rhoa.^2).*(1 - (2*(y.^2 + a^2)./rhoa.^2));
ISigmaXY =   prim * (x./rhoa.^2).*(1 - ((2*(y.^2))./rhoa.^2));
ISigmaZZ =  -prim * (2*psi) .* (y./rhoa.^2).*(1 + (a^2./rhoa.^2));
%%%SigmaZZ = psi*(SigmaXX + SigmaYY)

%%%%%%     Self Stress  Screw (Cai {Non-Singular}
% SigmaXZ = -(miu*b)/(2*pi) * (y./rhoa.^2) .* (1 + (a^2/rhoa.^2));
% SigmaYZ = -(miu*b)/(2*pi) * (x./rhoa.^2) .* (1 + (a^2/rhoa.^2));
% SELF STRESS PLOT
x=x+l;  % acomodate the plot x axe
h2=figure(5);
subplot(2,2,1); 
contourf(x, y, ISigmaZZ,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Im ZZ (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');
colorbar;
%caxis([-10 10])

subplot(2,2,2); 
contourf(x, y, ISigmaXX,50, 'ShowText', 'off', 'LineColor', 'none')
title('Im XX (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold');
colorbar;
%caxis([-15 15])

subplot(2,2,3); 
contourf(x, y, ISigmaYY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Im YY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold'); 
%caxis([-15 15])
colorbar

subplot(2,2,4); 
contourf(x, y, ISigmaXY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('Im XY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal
set(gca,'FontWeight','bold'); 
colorbar;
%caxis([-15 15])

set(h2,'position',[10 10 1000 1000]); colormap(jet);





function MakeMaps3 (S, Xm, Ym, G, param3)
nombres = ['Sxx'; 'Syy' ;'Szz'; 'Syz'];
maping = {'Airy', 'Self_Stress', 'Weygand_Stress','FEM + Self_Stress+ Weygang', 'SelfStress + ImageStress',...
        'Airy + SelfStress + ImageStress'};

    function plotmap1 % this function specify the caxis  :)
        if param3(1) == 0

            else
                caxis([param3(1) param3(2)])
        end
    end

    for o=1:4
        nombres(o,:);
        fig = figure;
        fig.Units  = 'centimeters';
        fig.Position(3) = 6;
        fig.Position(4) = 6;
        set(fig.Children,'FontName','Times','FontSize',10);
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space
        fig.PaperPositionMode   = 'auto';
        
        contourf(Xm,Ym, S(:,:,o), 150, 'ShowText', 'off', 'LineColor', 'none')
        %title(['',num2str(maping{G}),'',num2str(nombres(o,:)),'']); 
        colorbar
        %xlabel('X (µm)'); ylabel('Y (µm)'); 
        axis image
        set(gca,'FontWeight','bold');
        plotmap1
        colormap(jet)
       
        num2str(maping{G})
        print('-dpng','-r1000', ['',num2str(maping{G}),'-',num2str(nombres(o,:)),''])
    end
end