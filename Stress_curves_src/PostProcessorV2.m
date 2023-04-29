% This fucntion plots the internal stress and the total stress coming from paraview
% file extracted along a line over a free surface. several protting options
% Do not include the theorethical aspects.

clc
close all
clear all


A = importdata('/home/gonzalez/Desktop/BENCHMARKS/PaperBenchmarks/InfEdgeImproved/X120Y240Z20.csv');
%plotforpaper()
tiempo = A.data(:,3)/10000;

%%%%  FEM Data --------------------------------------------------
% Load of the FEM stress 
FS = [A.data(:,23) A.data(:,24) A.data(:,25) A.data(:,26) A.data(:,27) A.data(:,28)];
%FS = fempar2HandLothe(FS); %(Oxx, Oyy, Ozz, Oxy, Oxz, Oyz)
FS=FS/1e6;
%MakeMaps1(FS(:,1), FS(:,2), FS(:,3), FS(:,4), tiempo, 2,1);
%MakeMapsAll(FS(:,1), FS(:,2), FS(:,3), FS(:,4), FS(:,5), FS(:,6), tiempo, 1,1);


% Load of the Internal stress (Self or Weygand)
IS = [A.data(:,11) A.data(:,12) A.data(:,13) A.data(:,14) A.data(:,15) A.data(:,16)];
%IS = fempar2HandLothe(IS); %(Oxx, Oyy, Ozz, Oxy, Oxz, Oyz)
IS=IS/1e6;
%MakeMaps1(IS(:,1), IS(:,2), IS(:,3), IS(:,6), tiempo, 1,2);
%MakeMapsAll(IS(:,1), IS(:,2), IS(:,3), IS(:,4), IS(:,5), IS(:,6), tiempo, 1,2);

% suma
%MakeMapsAll(FS(:,1)+IS(:,1), FS(:,2)+IS(:,2), FS(:,3)+IS(:,3), FS(:,4)+IS(:,4), FS(:,5)+IS(:,5), FS(:,6)+IS(:,6), tiempo, 1,3);


[psi, b, miu, Y, X, l] = material;  % Material properties for analytics

%%%% Theorethical part------------------------------------------
% Airy
[Aixx,Aiyy,Aixy,Aizz] = AiryStress(b, miu, l, psi, X, tiempo);
%MakeMaps1(Aixx,Aiyy,Aixy,Aizz, tiempo,2,1)

% Image
[ImSxx,ImSyy,ImSxy,ImSzz] = SelfStressCai(psi, -b, miu, X-l, tiempo);
%MakeMaps1(ImSxx,ImSyy,ImSxy,ImSzz, tiempo,1,2)

% Selfstress
[Sxx,Syy,Sxy,Szz] = SelfStressCai(psi,b,miu,X+l, tiempo);
%MakeMaps1(Sxx,Syy,Sxy,Szz,tiempo,1,1)

% airy + Image
%MakeMaps1(ImSxx+Aixx,ImSyy+Aiyy,ImSxy+Aixy,ImSzz+Aizz,tiempo,1,1)

%savefigure(list(i).name, 1)
%close all 



%----------------------- FIGURES PAPER-------------------------------------
plotforpaper()
% Theorethical Self Stress 
MakeMaps2(Sxx,Sxy, tiempo,1,1) 

%El-numodis Selfstress
ps = 2
MakeMaps2(IS(1:ps:end,1),-IS(1:ps:end,2), tiempo(1:ps:end), 2,1); 

% El-Numodis Total Stress
ps=12
MakeMaps2(FS(1:ps:end,1)+IS(1:ps:end,1), FS(1:ps:end,2)+IS(1:ps:end,2), tiempo(1:ps:end), 2,1);
% theoretical total stress
MakeMaps2(ImSxx(1:ps:end)+Aixx(1:ps:end)+Sxx(1:ps:end),ImSxy(1:ps:end)+Aixy(1:ps:end)+Sxy(1:ps:end),tiempo(1:ps:end),4,1)
xlim([-1 1])
legend('H&L \sigmaxx','H&L \sigmaxy','Sim \sigmaxx','Sim \sigmaxy', 'location','southwest')
 

disp("Finished")
%%%%%%%%%%%%  FUNCTIONS  %%%%%%%%%%%%%%%

function MakeMaps2(SigmaXX,SigmaXY, y,G,h)
figure(h) 
    if G==1
    plot(y,SigmaXX,'k-','linewidth', 1)
    hold on
    plot(y,SigmaXY,'r-','linewidth', 1)%, 'linewidth', 1.5)
%     legend('Self \sigma_{XX}','Self \sigma_{XY}')
    end
    if G==2
    hold on
    plot(y,SigmaXX,'ko', 'markersize',2,'linewidth', 1)
    plot(y,SigmaXY,'ro','markersize',2,'linewidth', 1)
    %legend('\sigma_{XX}','\sigma_{XY}')
    end
    
    if G==3
    hold on
    plot(y,SigmaXX,'k--','linewidth', 1)
    plot(y,SigmaXY,'r--','linewidth', 1)
    %legend('\sigma_{XX}','\sigma_{XY}')
    end
    
    if G==4
    hold on
    plot(y,SigmaXX,'k-','linewidth',1)
    plot(y,SigmaXY,'r-','linewidth',1)
%     legend('Total \sigma_{XX}','Total \sigma_{XY}')
    end
    
    %hold off
    legend box off
    grid on
    set(gca,'FontWeight','bold');
    set(gca,'linewidth',1)
    xlabel('X (\mum)'); ylabel('\sigma (MPa)'); 
    %xlim([-0.5 0.5]);
end


%%%-----------------------------------------------------------------------
function plotforpaper()
    fig = figure;
    fig.Units  = 'centimeters';
    fig.Position(3) = 9;
    fig.Position(4) = 9;
    set(fig.Children,'FontName','Times','FontSize',10);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space
    fig.PaperPositionMode   = 'auto';
    axis square
    % %set(gcf, 'color', 'none'); set(gca, 'color', 'none');
end %----------------------------------------------------------------------


function MakeMaps1(SigmaXX,SigmaYY,SigmaXY,SigmaZZ, y,G,h)
figure(h)
 
    if G==1
    plot(y,SigmaXX,'k-','linewidth', 2)
    hold on
    plot(y,SigmaYY,'b-','linewidth', 2)
    plot(y,SigmaXY,'r-','linewidth', 2)
    plot(y,SigmaZZ,'g-','linewidth', 2)
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}')
    end
    if G==2
    hold on
    plot(y,SigmaXX,'.k')
    plot(y,SigmaYY,'b.')
    plot(y,SigmaXY,'r.')
    plot(y,SigmaZZ,'g.')
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}','Theory')
    end
    
    %hold off
    legend box off
    grid on
    set(gca,'FontWeight','bold');
    set(gca,'linewidth',1)
    xlabel('X (\mum)'); %ylabel('\sigma (MPa)'); 
    %xlim([-1 1]);
    
%      ylim([-4 4])

    %print('-dpdf','-r600', 'AlongLineBoundary')
    %print('-djpeg','-r600', 'statJon_BoundLineRefWeygand')  
end
%%%%-----------------------------------------------------------------------

% this function convert the components coming from the paraview to H&L
function [Oxx, Oyy, Oxy, Ozz, Oxz, Oyz] = matchHL(xx, yy, zz, xy, xz, yz)
    Oxx = yy;  Oyy = zz;  Oxy = yz;  Ozz = xx;  Oxz = xy;  Oyz = xz;
end

function [ostr] = fempar2HandLothe(istr)
    ostr(:,1) = istr(:,2); ostr(:,2) = istr(:,3); ostr(:,3) = istr(:,6);
    ostr(:,4) = istr(:,1); ostr(:,5) = istr(:,4); ostr(:,6) = istr(:,5);
end


function MakeMapsAll(SigmaXX,SigmaXY,SigmaXZ,SigmaYY,SigmaYZ,SigmaZZ, y,G,h)
figure(h)

    if G==1
    plot(y,SigmaXX,'k-','linewidth', 2)
    hold on
    plot(y,SigmaXY,'b-','linewidth', 2)
    plot(y,SigmaXZ,'c-','linewidth', 2)
    plot(y,SigmaYY,'g-','linewidth', 2)
    plot(y,SigmaYZ,'y-','linewidth', 2)
    plot(y,SigmaZZ,'r-','linewidth', 2)
    legend('\sigma_{XX}', '\sigma_{XY}','\sigma_{XZ}','\sigma_{YY}','\sigma_{YZ}','\sigma_{ZZ}')
    end
    %hold off
    legend box off
    grid on
    set(gca,'FontWeight','bold');
    set(gca,'linewidth',1)
    xlabel('X (\mum)'); %ylabel('\sigma (MPa)'); 
    %xlim([-1 1]);   
%      ylim([-4 4])

    %print('-dpdf','-r600', 'AlongLineBoundary')
    %print('-djpeg','-r600', 'statJon_BoundLineRefWeygand')  
end
%%%%-----------------------------------------------------------------------


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
end

% parameter for gererate theoretical fields
function [psi, b, miu, Y, X, l] = material()
% generate the theorethical Airy with the absice of imput_data
% % Material properties: Copper (Numodis)
    %E = 111216;            % Young Modulus  [Mpa]
    psi = 0.324;             % Poisson ratio
    b = 0.2552;              % Burgers vector [nm]
    miu = 42000;             % [MPa]
    %E = 2*miu*(1+psi);
    b = b * 0.001;     % nm    ==>  micrometer

    Y = [-1:0.01:1];  % micrometers
    X = 0;
    l =0.1;
end  %---------------------------------------------------------------------

