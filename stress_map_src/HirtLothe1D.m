clc
clear all
close all

% % Material properties: Copper (Numodis)
%E = 111216 ;            % Young Modulus  [Mpa]
psi = 0.324;             % Poisson ratio
b = 0.2552;              % Burgers vector [nm]
miu = 42000;             % [MPa]
%E = 2*miu*(1+psi);
b = b * 0.001;     % nm    ==>  micrometer

Y = [-1:0.01:1];  % micrometers
X = 0;
l =0.1; 


fig = figure;
fig.Units  = 'centimeters';
fig.Position(3) = 9;
fig.Position(4) = 9;
set(fig.Children,'FontName','Times','FontSize',10);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space
fig.PaperPositionMode   = 'auto';
axis square
%set(gcf, 'color', 'none'); set(gca, 'color', 'none');



% Selfstress
[Sxx,Syy,Sxy,Szz] = SelfStressCai(psi,b,miu,X+l,Y);
%MakeMaps1(Sxx,Syy,Sxy,Szz,Y,1)
%[Sxx,Syy,Sxy,Szz] = SelfStressHirt(psi,b,miu,X+l,Y);
%MakeMaps1(Sxx,Syy,Sxy,Szz,Y,2)

% Image_Selfstress
[ImSxx,ImSyy,ImSxy,ImSzz] = SelfStressCai(psi,-b,miu,X-l,Y);
%MakeMaps1(ImSxx,ImSyy,ImSxy,ImSzz, Y,1)

% Airy function
[Aixx,Aiyy,Aixy,Aizz] = AiryStress(b,miu,l,psi,X,Y);
%MakeMaps1(Aixx,Aiyy,Aixy,Aizz, Y,1)

% Image_Selfstress + Airy
%MakeMaps1(Aixx+ImSxx,Aiyy+ImSyy,Aixy+ImSxy,Aizz+ImSzz, Y,1)

% Image_Selfstress + Selfstress
%MakeMaps1(Sxx+ImSxx,Syy+ImSyy,Sxy+ImSxy,Szz+ImSzz, Y,1)

% Selfstress + Airy
%MakeMaps1(Aixx+Sxx,Aiyy+Syy,Aixy+Sxy,Aizz+Szz, Y,1)

% Total stress =  Self + Image + Airy
MakeMaps1(Sxx+ImSxx+Aixx,Syy+ImSyy+Aiyy,Sxy+ImSxy+Aixy,Szz+ImSzz+Aizz, Y,1)








%A = importdata('../SelfStressInf2/FEMBoundaryRefinement.csv');
%MakeMaps1(A.data(:,20)*1e-6,A.data(:,22)*1e-6,A.data(:,21)*1e-6,A.data(:,17)*1e-6, A.data(:,4)/10000,2)


disp FINISHED

%%%%-----------------------------------------------------------------------
%%%%%%      Self Stress  Edge (Hirt and Lothe {Singular})          %%%%%%%%
%%%%-----------------------------------------------------------------------
function [SigmaXX,SigmaYY,SigmaXY,SigmaZZ]=SelfStressHirt(psi,b,miu,x,y)
prim = (miu*b)/(2*pi*(1-psi));
SigmaXX =  -prim * ( (y.*(3*x.^2 + y.^2))./ (x.^2 + y.^2).^2 );
SigmaYY =   prim * ( (y.*(x.^2 - y.^2))./ (x.^2 + y.^2).^2 );
SigmaXY =   prim * ( (x.*(x.^2 - y.^2))./ (x.^2 + y.^2).^2 );
SigmaZZ =  -prim *(2*psi) * (y./(x.^2 + y.^2)); 
%SigmaZZ = psi*(SigmaXX + SigmaYY);
%%%%%%     Self Stress  Screw (Hirt and Lothe {Singular})
% SigmaXZ = -(miu*b)/(2*pi) * (y./(x.^2 + y.^2));
% SigmaYZ = -(miu*b)/(2*pi) * (x./(x.^2 + y.^2));
end

%%%%-----------------------------------------------------------------------
%%%%%%           Self Stress  Edge (Cai {Non-Singular})            %%%%%%%%
%%%%-----------------------------------------------------------------------
function [SigmaXX,SigmaYY,SigmaXY,SigmaZZ]=SelfStressCai(psi,b,miu,x,y)
a = 0.0003;
prim = (miu*b)/(2*pi*(1-psi));
rhoa = sqrt(a^2 + x^2 + y.^2);
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
%%%%%%      Airy Stress Edge (Hirt and Lothe {Singular}) 
%%%%-----------------------------------------------------------------------
function [I_SigmaXX,I_SigmaYY,I_SigmaXY,I_SigmaZZ]=AiryStress(b,miu,l,psi,x,y)
l = l;
r = sqrt(l.^2 + y.^2);
prim = (miu*b)/(2*pi*(1-psi));
I_SigmaXX = -prim * (4*l*x.*y)./(r.^6) .* (3.*(l-x).^2 - y.^2);
I_SigmaYY =  prim * (2*l)./(r.^6) .* ( 4*y.*(l-x).^3 + 6*x.*y.*(l-x).^2 + 4*(y.^3).*(l-x) - 2*x.*y.^3);
I_SigmaXY = -prim * (2*l)./(r.^6) .* ( (l-x).^4 + 2*x.*(l-x).^3 - 6*x.*(y.^2).*(l-x) - y.^4);
I_SigmaZZ =  prim * (8*l*psi)./(r.^6) .* (y.*(l-x).^3 + (l-x).*(y).^3); 
%%%%-----------------------------------------------------------------------
end


function MakeMaps1(SigmaXX,SigmaYY,SigmaXY,SigmaZZ, y,G)

    
    if G==1
    plot(y,SigmaXX,'k-','linewidth', 2)
    hold on
    plot(y,SigmaYY,'b-','linewidth', 2)
    plot(y,SigmaXY,'r-','linewidth', 2)
    plot(y,SigmaZZ,'g-','linewidth', 2)
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}')
    end
    if G==2
    plot(y,SigmaXX,'.k')
    plot(y,SigmaYY,'b.')
    plot(y,SigmaXY,'r.')
    plot(y,SigmaZZ,'g.')
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}','FEM')
    end
    
    
    %hold off
    legend box off
    grid on
    set(gca,'FontWeight','bold');
    set(gca,'linewidth',1)
    xlabel('Y (\mum)'); ylabel('\sigma (MPa)'); 

    xlim([-1 1])
    %print('-dpdf','-r600', 'AlongLineBoundary')
    %print('-djpeg','-r600', 'LineBoundaryAiry')



  
end