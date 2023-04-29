%%%------------------------------------------------------------------------
%  function to plot the Self Stress, Image Stress, Irish function and 
%  the Sum of some of the previous parameters.
%
% Function usage:
%
%   MakeMaps1(Sxx,Syy,Sxy,Szz, X_Meshgrid, Y_meshgrid, FigureNumber, 
%             Colorbarlimits[-10,10])
%   -- If colorbarnumber is [0 0] the colorbar takes free range :)       
%
%%%------------------------------------------------------------------------


clc
clear all
close all

% % Material properties: Copper (Numodis)
%E = 111216 ;            % Young Modulus  [Mpa]
% psi = 0.324;             % Poisson ratio
% b = 0.2552;              % Burgers vector [nm]
% miu = 42000;             % [MPa]
% E = 2*miu*(1+psi);
% b = b * 0.001;     % nm    ==>  micrometer

l = 1e-06
miu = 42e9
b = 2.552e-10
psi = 0.324


Y = [-4:0.025:4]*1e-6;  % micrometers
X = [-4:0.025:0]*1e-6;


% Selfstress#
X1 = X+l;  % micrometers
[x,y]=meshgrid(X1,Y);
[Sxx,Syy,Sxy,Szz] = SelfStressCai(psi,b,miu,x,y);
MakeMaps1(Sxx,Syy,Sxy,Szz, x-1, y, 1, [0 0])


% Image_Selfstress
X1 = X-l;  % micrometers
[x,y]=meshgrid(X1,Y);
[ImSxx,ImSyy,ImSxy,ImSzz] = SelfStressCai(psi,-b,miu,x,y);
MakeMaps1(ImSxx,ImSyy,ImSxy,ImSzz, x+1, y, 2, [0,0])

% Airy function
[x,y]=meshgrid(X,Y);
[Aixx,Aiyy,Aixy,Aizz] = AiryStress(b,miu,l,psi,x,y);
MakeMaps1(Aixx,Aiyy,Aixy,Aizz, x, y, 3, [0,0])

% Image_Selfstress + SelfStress
[x,y]=meshgrid(X,Y);
MakeMaps1(ImSxx+Sxx,ImSyy+Syy,ImSxy+Sxy,ImSzz+Szz, x, y, 4, [0 0])

% Image_Selfstress + SelfStress + Airy
[x,y]=meshgrid(X,Y);
MakeMaps1(Aixx+ImSxx+Sxx,Aiyy+ImSyy+Syy,Aixy+ImSxy+Sxy,Aizz+ImSzz+Szz, x, y, 5, [0 0])

% Image_Selfstress + Airy
[x,y]=meshgrid(X,Y);
MakeMaps1(Aixx+ImSxx,Aiyy+ImSyy,Aixy+ImSxy,Aizz+ImSzz, x, y, 6, [0, 0])

%%%%-----------------------------------------------------------------------
%%%%%%           Self Stress  Edge (Cai {Non-Singular})            %%%%%%%%
%%%%-----------------------------------------------------------------------
function [SigmaXX,SigmaYY,SigmaXY,SigmaZZ]=SelfStressCai(psi,b,miu,x,y)
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
end



%%%%-----------------------------------------------------------------------
%%%%%%      Total Stress Edge (Hirt and Lothe {Singular}) 
%%%%-----------------------------------------------------------------------
function [I_SigmaXX,I_SigmaYY,I_SigmaXY,I_SigmaZZ]=AiryStress(b,miu,l,psi,x,y)
l = l;
r = sqrt((x-l).^2 + y.^2);
prim = (miu*b)/(2*pi*(1-psi));
I_SigmaXX = -prim * (4*l*x.*y)./(r.^6) .* (3.*(l-x).^2 - y.^2);
I_SigmaYY =  prim * (2*l)./(r.^6) .* ( 4*y.*(l-x).^3 + 6*x.*y.*(l-x).^2 + 4*(y.^3).*(l-x) - 2*x.*y.^3);
I_SigmaXY = -prim * (2*l)./(r.^6) .* ( (l-x).^4 + 2*x.*(l-x).^3 - 6*x.*(y.^2).*(l-x) - y.^4);
I_SigmaZZ =  prim * (8*l*psi)./(r.^6) .* (y.*(l-x).^3 + (l-x).*(y).^3); 
%%%%-----------------------------------------------------------------------
end



function MakeMaps1(SigmaXX,SigmaYY,SigmaXY,SigmaZZ, x, y, G, param3)


    function plotmap1 % this function specify the caxis  :)
        if param3(1) == 0

        else
            caxis([param3(1) param3(2)])
        end
    end

    h1=figure(G);
    subplot(2,2,1); 
    imagesc(x(1,:), y(:,1), SigmaZZ); set(gca,'YDir','normal')
    %contourf(x, y, SigmaZZ,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
    title('ZZ (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');colorbar;
    plotmap1

    subplot(2,2,2); 
    imagesc(x(1,:), y(:,1),SigmaXX); set(gca,'YDir','normal')
    %contourf(x, y, SigmaXX,50, 'ShowText', 'off', 'LineColor', 'none')
    title('XX (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');colorbar;
    plotmap1

    subplot(2,2,3); 
    imagesc(x(1,:), y(:,1), SigmaYY); set(gca,'YDir','normal')
    %contourf(x, y, SigmaYY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
    title('YY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');
    colorbar
    plotmap1

    subplot(2,2,4); 
    imagesc(x(1,:), y(:,1), SigmaXY); set(gca,'YDir','normal')
    %contourf(x, y, SigmaXY,50, 'ShowText', 'off', 'LineColor', 'none')%;hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
    title('XY (Mpa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold'); colorbar;
    
    plotmap1

    set(h1,'position',[10 10 1000 1000]); colormap(jet);
    
    names = {'SelfStress', 'ImageStress','Airy function', 'SelfStress + ImageStress',...
            'Airy + SelfStress + ImageStress', 'Airy+ ImageStress'};
    suptitle(names(G));
    
%    if G == 5      
%    print('-dpng','-r600', G)
%    end
end