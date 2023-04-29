%%%%
%%%%
%%%%

clc
close all
clear all

[psi, b, miu, Y, X, l] = material;  % Material properties for analytics


%%%%%%   list the *.csv files on thhe folder
% cd ../SelfStressInf2/;
list = dir('Si_Weygand/*StressPar*.csv')
length(list)

for i = 1:length(list)
    A = importdata(['No_Weygand/',num2str(list(i).name),'']);
    %plotforpaper()
    tiempo = A.data(:,28)/10000;
    
    %%%%  FEM Data --------------------------------------------------
    % Load of the FEM stress 
    FS = [A.data(:,1) A.data(:,2) A.data(:,3) A.data(:,4) A.data(:,5) A.data(:,6)];
    FS = fempar2HandLothe(FS); %(Oxx, Oyy, Ozz, Oxy, Oxz, Oyz)
    FS=FS/1e6;
    MakeMaps1(FS(:,1), FS(:,2), FS(:,6), FS(:,4), tiempo, 2,1);
    
    % Load of the Internal stress (Self or Weygand)
    IS = [A.data(:,14) A.data(:,15) A.data(:,16) A.data(:,17) A.data(:,18) A.data(:,19)];
    IS = fempar2HandLothe(IS); %(Oxx, Oyy, Ozz, Oxy, Oxz, Oyz)
    IS=IS/1e6;
    MakeMaps1(IS(:,1), IS(:,2), IS(:,3), IS(:,4), tiempo, 1,2);
    
    % FEM+ (Self or Weygand)
    MakeMaps1(FS(:,1)+IS(:,1), FS(:,2)+IS(:,2), FS(:,6)+IS(:,3), FS(:,4)+IS(:,4), tiempo, 1,3);
    
    
    %%%% Theorethical part------------------------------------------
    % Airy
    [Aixx,Aiyy,Aixy,Aizz] = AiryStress(b, miu, l, psi, X, tiempo);
    MakeMaps1(Aixx,Aiyy,Aixy,Aizz, tiempo,2,7)

    % Image
    [ImSxx,ImSyy,ImSxy,ImSzz] = SelfStressCai(psi, -b, miu, X-l, tiempo);
    MakeMaps1(ImSxx,ImSyy,ImSxy,ImSzz, tiempo,1,2)

    % Selfstress
    [Sxx,Syy,Sxy,Szz] = SelfStressCai(psi,b,miu,X+l, tiempo);
    MakeMaps1(Sxx,Syy,Sxy,Szz,tiempo,1,1)

    % airy + Image
    MakeMaps1(ImSxx+Aixx,ImSyy+Aiyy,ImSxy+Aixy,ImSzz+Aizz,tiempo,1,1)
    
    % airy + self
    MakeMaps1(Sxx+Aixx,Syy+Aiyy,Sxy+Aixy,Szz+Aizz,tiempo,1,2)
    
    % Image +self
    MakeMaps1(ImSxx+Sxx,ImSyy+Syy,ImSxy+Sxy,ImSzz+Szz,tiempo,1,4)

    % airy + Image +self
    MakeMaps1(ImSxx+Aixx+Sxx,ImSyy+Aiyy+Syy,ImSxy+Aixy+Sxy,ImSzz+Aizz+Szz,tiempo,1,5)


    %savefigure(list(i).name, 1)
    %close all
    %%%% Statistical part----------------------------------------
    %plotforpaper()

     MakeMaps1(Djon(Aixx,FS(:,1)), Djon(Aiyy,FS(:,2)), ...
              Djon(Aixy,FS(:,6)), Djon(Aizz,FS(:,4)), tiempo,2, 1)
    ylim([-4 4])
%     
%     MakeMaps1(Djon(Aixx,FS(:,1)+IS(:,1)), Djon(Aiyy,FS(:,2)+IS(:,2)), ...
%               Djon(Aixy,FS(:,6)-IS(:,3)), Djon(Aizz,FS(:,4)+IS(:,4)), tiempo,2, 1)
%     ylim([-4 4])
    
    %savefigure(list(i).name, 2)

    
    %legend([h1 h2 h3 h5],{'h1','h2','h3','h5})
    %hLeg = legend('example')
    %set(hLeg,'visible','off')
    
 close all
end

disp FINISHED








%%%%%%%%%%%%  FUNCTIONS  %%%%%%%%%%%%%%%


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
    xlim([-1 1]);
    
%      ylim([-4 4])

    %print('-dpdf','-r600', 'AlongLineBoundary')
    %print('-djpeg','-r600', 'statJon_BoundLineRefWeygand')  
end
%%%%-----------------------------------------------------------------------

function out_Djon = Djon(Theo, Exp)
    if length(Theo) == length(Exp)
        out_Djon = ((Theo-Exp)./Theo);
    else
        disp('Observed and predicted values do not have same size...')
        disp('Aborting D_Jon function... ')
    end
end%-----------------------------------------------------------------------

%%%-----------------------------------------------------------------------
% reshape the six component stress tensor (in column format) to a
% meshgrided format to be plotted. 2 ways can be used Numodis format and
% ElmerSpreadsheet format
%%%-----------------------------------------------------------------------
function [S1,S2,S3,S4,S5,S6] = SplitStressComp(tobesplit, modo, disc1, disc2)
    switch modo
        case 1 % numodis format
            S1 = tobesplit(:,1); 
            S2 = tobesplit(:,2);
            S3 = tobesplit(:,3);
            S4 = tobesplit(:,4);
            S5 = tobesplit(:,5);
            S6 = tobesplit(:,6);
        case 2
            S1 = reshape(tobesplit(:,1), [disc2,disc1]);
            S2 = reshape(tobesplit(:,4), [disc2,disc1]);
            S3 = reshape(tobesplit(:,6), [disc2,disc1]);
            S4 = reshape(tobesplit(:,2), [disc2,disc1]);
            S5 = reshape(tobesplit(:,3), [disc2,disc1]);
            S6 = reshape(tobesplit(:,5), [disc2,disc1]);
    end
end


function savetxt(pos)
    disp '*** Saving point coordinate data for numodis Couplig...';
    disp '*** Saving mode: total#  x1 y1 z1 x2 y2 z2 x3 y3 z3 ...';
    point4numodis = reshape(pos(:,1:3)', [1 3*length(pos(:,1))]);
    fid1=fopen('pos.txt','w');
    fprintf(fid1, '%d \n', length(point4numodis));
    fprintf(fid1, '%f ', point4numodis');
    fclose(fid1);

end


% this function convert the components coming from the paraview to H&L
function [Oxx, Oyy, Oxy, Ozz, Oxz, Oyz] = matchHL(xx, yy, zz, xy, xz, yz)
    Oxx = yy;  Oyy = zz;  Oxy = yz;  Ozz = xx;  Oxz = xy;  Oyz = xz;
end

function [ostr] = fempar2HandLothe(istr)
    ostr(:,1) = istr(:,2); ostr(:,2) = istr(:,3); ostr(:,3) = istr(:,6);
    ostr(:,4) = istr(:,1); ostr(:,5) = istr(:,4); ostr(:,6) = istr(:,5);
end

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


function savefigure(este,n)

    nombre = este(1:length(este)-4)
    
    if n == 1
       fig_name = strcat('WAiry',num2str(nombre))

    else
       fig_name = strcat('WStat',num2str(nombre));

    end

     print(fig_name,'-djpeg','-r1000'); 
end


function MakeMaps2(SigmaXX,SigmaXY, y,G,h)
figure(h)
 
    if G==1
    plot(y,SigmaXX,'k-')
    hold on
    plot(y,SigmaXY,'r-')%, 'linewidth', 1.5)
    %legend('W \sigma_{XX}','W \sigma_{XY}')
    end
    if G==2
    hold on
    plot(y,SigmaXX,'k.')
    plot(y,SigmaXY,'b.')
    %legend('\sigma_{XX}','\sigma_{XY}')
    end
    
    %hold off
    legend box off
    grid on
    set(gca,'FontWeight','bold');
    set(gca,'linewidth',1)
    xlabel('X (\mum)'); ylabel('\sigma (MPa)'); 
    %xlim([-0.5 0.5]);
end