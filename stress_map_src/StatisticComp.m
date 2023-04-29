clc 
close all 
clear all


%%%%%%   list the *.csv files on thhe folder
cd ../SelfStressInf2/;
list = dir('*Line*.csv');
length(list);

[psi, b, miu, Y, X, l] = material;  % Material properties for analytics


for i = 1: length(list)    % looping for each eleements in the list
    
    Refi = importdata(['../SelfStressInf2/',num2str(list(i).name),'']);
    % FEM data
    Refi = [Refi.data(:,4) Refi.data(:,20) Refi.data(:,22) Refi.data(:,21) Refi.data(:,17)];
    % H&L stuffs       % Airy function
    [Aixx,Aiyy,Aixy,Aizz] = AiryStress(b, miu, l, psi, X, Refi(:,1)/10000);

    if contains( list(i).name,'Weygand')   %   Weigand case

        % MakeMaps1(Refi(:,2)/1e6, Refi(:,3)/1e6, Refi(:,4)/1e6, ...
        %           Refi(:,5)/1e6, Refi(:,1)/10000, 1,1);
        % MakeMaps1(Aixx,Aiyy,Aixy,Aizz, Refi(:,1)/10000,2,1)
        
        %%%% statistical part
       
        MakeMaps1(Djon(Aixx,Refi(:,2)/1e6), Djon(Aiyy,Refi(:,3)/1e6), ...
                  Djon(Aixy,Refi(:,4)/1e6), Djon(Aizz,Refi(:,5)/1e6), Refi(:,1)/10000,2, i)
        ylim([-4 4])
        title(list(i).name)
    
    else                                %   No Weigand case
        % FEM 
        % MakeMaps1(Refi(:,2)/1e6, Refi(:,3)/1e6, Refi(:,4)/1e6, ... 
        %            Refi(:,5)/1e6, Refi(:,1)/10000, 1,1);
        
        % Image_Selfstress
        [ImSxx,ImSyy,ImSxy,ImSzz] = SelfStressCai(psi, -b, miu, X-l, Refi(:,1)/10000);
        %MakeMaps1(ImSxx,ImSyy,ImSxy,ImSzz, Refi(:,1)/10000,1,2)
        
        MakeMaps1(Djon(Aixx+ImSxx,Refi(:,2)/1e6), Djon(Aiyy+ImSyy,Refi(:,3)/1e6), ...
                  Djon(Aixy+ImSxy,Refi(:,4)/1e6), Djon(Aizz+ImSzz,Refi(:,5)/1e6), Refi(:,1)/10000,2, i)
        ylim([-4 4])
        title(list(i).name)
    end
     
end





%plotforpapaer    % make firgure with publicable dimentions
%[psi, b, miu, Y, X, l] = material()   % Material properties

stop




% Selfstress
%[Sxx,Syy,Sxy,Szz] = SelfStressCai(psi,b,miu,X+l,Y);
%MakeMaps1(Sxx,Syy,Sxy,Szz,Y,1)
%[Sxx,Syy,Sxy,Szz] = SelfStressHirt(psi,b,miu,X+l,Y);
%MakeMaps1(Sxx,Syy,Sxy,Szz,Y,2)

% Image_Selfstress
%[ImSxx,ImSyy,ImSxy,ImSzz] = SelfStressCai(psi,-b,miu,X-l,Refi(:,1)/10000);
%MakeMaps1(ImSxx,ImSyy,ImSxy,ImSzz, Refi(:,1)/10000,1,2)

% Airy function
%[Aixx,Aiyy,Aixy,Aizz] = AiryStress(b,miu,l,psi,X,Refi(:,1)/10000);
%MakeMaps1(Aixx,Aiyy,Aixy,Aizz, Refi(:,1)/10000,2,1)

%MakeMaps1(Aixx+ImSxx,Aiyy+ImSyy,Aixy+ImSxy,Aizz+ImSzz, Refi(:,1)/10000,2,1)




% plot(Refi(:,1)/10000, Djon(Aiyy,Refi(:,3)/1e6))% D_Jon(Aixx,Refi(:,2)/1e6))
       
      
      
      
      
%%%%   FUNCTIONS  %%%%%%

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


function [psi, b, miu, Y, X, l] = material()
% generate the theorethical Airy with the absice of imput_data
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
end  %---------------------------------------------------------------------

function plotforpapaer()
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
    plot(y,SigmaXX,'k-')%,'linewidth', 2)
    hold on
    plot(y,SigmaYY,'b-')%,'linewidth', 2)
    plot(y,SigmaXY,'r-')%,'linewidth', 2)
    plot(y,SigmaZZ,'g-')%,'linewidth', 2)
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}')
    end
    if G==2
    hold on
    plot(y,SigmaXX,'.k')
    plot(y,SigmaYY,'b.')
    plot(y,SigmaXY,'r.')
    plot(y,SigmaZZ,'g.')
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}','Airy')
    end
    
    %hold off
    legend box off
    grid on
%     set(gca,'FontWeight','bold');
%     set(gca,'linewidth',1)
     xlabel('X (\mum)'); ylabel('\sigma (MPa)'); 
% 
%      xlim([-1 1])
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


function out_mae1 = mae1(obs, pred)
    if length(obs) == length(pred)
        N= length(obs)
        out_mae1 = 1/N *sum(abs(obs-pred))
    else
        disp('Observed and predicted values do not have same size...')
        disp('Aborting MAE function... ')
    end
end %-----------------------------------------------------------------------

function out_mse1 = mse1(Theo, Exp)
    if length(Theo) == length(Exp)
        N= length(Theo)
        out_mse1 = 1/N *sum((Theo-Exp)^2)
    else
        disp('Observed and predicted values do not have same size...')
        disp('Aborting MAE function... ')
    end
end%-----------------------------------------------------------------------

function out_rsquare1 = rsquare1(Theo, Exp)
    if length(Theo) == length(Exp)
        out_rsquare1 = 1 - (sum((Theo-Exp)^2)/(sum(Theo - mean(Exp))^2))
    else
        disp('Observed and predicted values do not have same size...')
        disp('Aborting MAE function... ')
    end
end%-----------------------------------------------------------------------

function plotstat1()
end