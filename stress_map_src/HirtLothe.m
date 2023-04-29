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


close all

% % Material properties: Copper (Numodis)
%E = 111216 ;            % Young Modulus  [Mpa]
psi = 0.324;             % Poisson ratio
b = 0.2552;              % Burgers vector [nm]
miu = 42000;             % [MPa]
E = 2*miu*(1+psi);
b = b * 0.001;     % nm    ==>  micrometer

Y = [-1:0.002:1];  % micrometers
X = [-1:0.002:0];

%Y = fictY/10000;
%X = fictX/10000;

l =0.1; 


% Selfstress
X1 = X+l;  % micrometers
[x,y]=meshgrid(X1,Y);
[Sxx,Syy,Sxy,Szz] = SelfStressCai(psi,b,miu,x,y);
%MakeMaps1(Sxx,Syy,Sxy,Szz, x-l, y, 1, [-20 20],l);
%MakeMaps3(cat(6,Sxx, Syy, Sxy, Szz), x, y, 5, [-20 20])


%[Sxx,Syy,Sxy,Szz] = SelfStressHirt(psi,b,miu,x,y);
%MakeMaps1(Sxx,Syy,Sxy,Szz, x-l, y, 2, [-20 20],l);
%MakeMaps3(cat(6,Sxx, Syy, Sxy, Szz), x, y, 4, [-20 20])



% Image_Selfstress
X1 = X-l;  % micrometers
[x,y]=meshgrid(X1,Y);
[ImSxx,ImSyy,ImSxy,ImSzz] = SelfStressCai(psi,-b,miu,x,y);
%MakeMaps1(ImSxx,ImSyy,ImSxy,ImSzz, x+l, y, 2, [-20,20],l)

% Airy function
[x,y]=meshgrid(X,Y);
[Aixx,Aiyy,Aixy,Aizz] = AiryStress(b,miu,l,psi,x,y);
MakeMaps1(Aixx,Aiyy,Aixy,Aizz, x, y, 3, [-20 20],l)

% Image_Selfstress + SelfStress
[x,y]=meshgrid(X,Y);
%MakeMaps1(ImSxx+Sxx,ImSyy+Syy,ImSxy+Sxy,ImSzz+Szz, x, y, 4, [-30 30],l)
%MakeMaps3(cat(6,ImSxx+Sxx,ImSyy+Syy,ImSxy+Sxy,ImSzz+Szz), x, y, 4, [-30 30])


% Image_Selfstress + SelfStress + Airy
[x,y]=meshgrid(X,Y);
%MakeMaps1(Aixx+ImSxx+Sxx,Aiyy+ImSyy+Syy,Aixy+ImSxy+Sxy,Aizz+ImSzz+Szz, x, y, 5, [0 0],l)

% Image_Selfstress + Airy
[x,y]=meshgrid(X,Y);
%MakeMaps1(Aixx+ImSxx,Aiyy+ImSyy,Aixy+ImSxy,Aizz+ImSzz, x, y, 6, [-20 20],l)

disp FINISHED




%%%-----------------------------------------------------------------------
function MakeMaps3 (S, Xm, Ym, G, param3)
nombres = ['Sxx'; 'Syy' ;'Szz'; 'Sxy'; 'Sxz'; 'Syz'];
maping = {'FemStress'; 'SelfStress'; 'FEM+SelfStress'; 'TheoHirtSelfStressEdge'; 'TheoCaiSelfStressEdge'};

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
        h = pcolor(S(:,:,o))
        set(h,'EdgeColor' ,'none')
        shading interp
        
        %contourf(Xm,Ym, S(:,:,o), 150, 'ShowText', 'off', 'LineColor', 'none')
        %title(['',num2str(maping{G,:}),'',num2str(nombres(o,:)),'']); 
        %colorbar
        %xlabel('X (µm)'); ylabel('Y (µm)'); 
        axis image
        axis off
        set(gca,'FontWeight','bold');
        plotmap1
        colormap(jet)
       
        %print('-dpng','-r1000', ['',num2str(maping{G,:}),'-',num2str(nombres(o,:)),''])
    end
end
%%%-----------


%%%%-----------------------------------------------------------------------
%%%%%%      Self Stress  Edge (Hirt and Lothe {Singular})          %%%%%%%%
%%%%-----------------------------------------------------------------------
function [SigmaXX,SigmaYY,SigmaXY,SigmaZZ]=SelfStressHirt(psi,b,miu,x,y)
prim = (miu*b)/(2*pi*(1-psi));
% % SigmaXX =  -prim * ( (y.*(3*x.^2 + y.^2))./ (x.^2 + y.^2).^2 );
% % SigmaYY =   prim * ( (y.*(x.^2 - y.^2))./ (x.^2 + y.^2).^2 );
% % SigmaXY =   prim * ( (x.*(x.^2 - y.^2))./ (x.^2 + y.^2).^2 );
% % SigmaZZ =  -prim *(2*psi) * (y./(x.^2 + y.^2)); 
%%SigmaZZ = psi*(SigmaXX + SigmaYY);
%%%%%%     Self Stress  Screw (Hirt and Lothe {Singular})
SigmaXX = -(miu*b)/(2*pi) * (y./(x.^2 + y.^2)); %Sxz
SigmaYY= -(miu*b)/(2*pi) * (x./(x.^2 + y.^2)); %Syz
SigmaXY = -(miu*b)/(2*pi) * (y./(x.^2 + y.^2)); %Sxz
SigmaZZ = -(miu*b)/(2*pi) * (x./(x.^2 + y.^2)); %Syz
end
%%%%-----------------------------------------------------------------------


%%%%-----------------------------------------------------------------------
%%%%%%           Self Stress  Edge (Cai {Non-Singular})            %%%%%%%%
%%%%-----------------------------------------------------------------------
function [SigmaXX,SigmaYY,SigmaXY,SigmaZZ]=SelfStressCai(psi,b,miu,x,y)
a = 0.000003;
prim = (miu*b)/(2*pi*(1-psi));
rhoa = sqrt(a^2 + x.^2 + y.^2);
SigmaXX =  -prim * (y./rhoa.^2).*(1 + (2*(x.^2 + a^2)./rhoa.^2));
SigmaYY =   prim * (y./rhoa.^2).*(1 - (2*(y.^2 + a^2)./rhoa.^2));
SigmaXY =   prim * (x./rhoa.^2).*(1 - ((2*(y.^2))./rhoa.^2));
SigmaZZ =  -prim * (2*psi) .* (y./rhoa.^2).*(1 + (a^2./rhoa.^2));
%%%SigmaZZ = psi*(SigmaXX + SigmaYY)
%%%%%%     Self Stress  Screw (Cai {Non-Singular}
%  SigmaXX = -(miu*b)/(2*pi) * (y./rhoa.^2) .* (1 + (a^2./rhoa.^2)); %Sxz
%  SigmaYY = -(miu*b)/(2*pi) * (x./rhoa.^2) .* (1 + (a^2./rhoa.^2)); %Syz
%  SigmaXY = -(miu*b)/(2*pi) * (y./rhoa.^2) .* (1 + (a^2./rhoa.^2)); %Sxz
%  SigmaZZ = -(miu*b)/(2*pi) * (x./rhoa.^2) .* (1 + (a^2./rhoa.^2)); %Syz
end
%%%%-----------------------------------------------------------------------



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
end
%%%%-----------------------------------------------------------------------


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
    imagesc(x(1,:), y(:,1), SigmaZZ); set(gca,'YDir','normal')
    %xline(0,'--','linewidth',2)
    hold on
    %plot(-l,0,'k.','markersize', 10)
    title('\sigma_{zz} (MPa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');colorbar;
    plotmap1

    subplot(2,2,2); 
    imagesc(x(1,:), y(:,1),SigmaXX); set(gca,'YDir','normal')
    %xline(0,'--','linewidth',2)
    hold on
    %plot(-l,0,'k.','markersize', 10)
    title('\sigma_{xx} (MPa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');colorbar;
    plotmap1

    subplot(2,2,3); 
    imagesc(x(1,:), y(:,1), SigmaYY); set(gca,'YDir','normal')
    %xline(0,'--','linewidth',2)
    hold on
    %plot(-l,0,'k.','markersize', 10)
    title('\sigma_{yy} (MPa)'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
    axis image
    set(gca,'FontWeight','bold');
    colorbar
    plotmap1

    subplot(2,2,4); 
    imagesc(x(1,:), y(:,1), SigmaXY); set(gca,'YDir','normal')
    %xline(0,'--','linewidth',2)
    hold on
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



