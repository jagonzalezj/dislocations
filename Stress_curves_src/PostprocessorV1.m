%%%%
%%%%
%%%%

clc
close all
clear all

[psi, b, miu, Y, X, l] = material;  % Material properties for analytics


%%%%%%   list the *.csv files on thhe folder
% cd ../SelfStressInf2/;
list = dir('Si_Weygand/*StressPar*.csv');
length(list);

for i = 1:length(list)
    %plotforpaper();
    list(i).name
    for tipo =1:2
        if tipo == 1
            % Weygand
            A = importdata(['/home/gonzalez/Desktop/BENCHMARKS/WEYGAND/Weygand_study2NonPBC/Si_Weygand/',num2str(list(i).name),'']);
        else
            % No Weygand
            A = importdata(['/home/gonzalez/Desktop/BENCHMARKS/WEYGAND/Weygand_study2NonPBC/No_Weygand/',num2str(list(i).name),'']);
        end
        
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
        MakeMaps1(IS(:,1), IS(:,2), IS(:,3), IS(:,4), tiempo, 2,1);
        
        % FEM-Image
        MakeMaps1(FS(:,1)+IS(:,1), FS(:,2)+IS(:,2), FS(:,6)-IS(:,3), FS(:,4)+IS(:,4), tiempo, 1,1);
        
        
        %%%% Theorethical part------------------------------------------
        % Airy
        [Aixx,Aiyy,Aixy,Aizz] = AiryStress(b, miu, l, psi, X, tiempo);
        MakeMaps1(Aixx,Aiyy,Aixy,Aizz, tiempo,2,1)
        
        % Image
        [ImSxx,ImSyy,ImSxy,ImSzz] = SelfStressCai(psi, -b, miu, X-l, tiempo);
        MakeMaps1(ImSxx,ImSyy,ImSxy,ImSzz, tiempo,1,2)
        
        % Selfstress
        [Sxx,Syy,Sxy,Szz] = SelfStressCai(psi,b,miu,X+l, tiempo);
        MakeMaps1(Sxx,Syy,Sxy,Szz,tiempo,1,1)
        
        % airy + Image
        %MakeMaps1(ImSxx+Aixx,ImSyy+Aiyy,ImSxy+Aixy,ImSzz+Aizz,tiempo,1,1)
        
        %savefigure(list(i).name, 1)
        %close all
        %%%% Statistical part----------------------------------------

        
        
        
        
        % version all togeteher 2222
        if i == 1 && tipo == 1
            plotforpaper()

            plot(tiempo, Aixy, 'r-','linewidth', 1)
            xlabel('X (\mum)'); ylabel('\sigma (MPa)'); 
            set(gca,'FontWeight','bold');
            set(gca,'linewidth',1)

            xlim([-0.2 0.2])
            hold on
        end
        simbolo = {'b--','b:', 'b-', 'k--','k:', 'k-'};
             
        if tipo == 1 % Weygand
            
            % version all togeteher 2222
            ps=1;
            %plot(tiempo(1:ps:end), FS(1:ps:end,6), simbolo{i},'linewidth', 1.5);
            plot(tiempo, FS(:,6), simbolo{i},'linewidth', 1);

            
            % version of the three graph separated  1111
            %figure(i)
            %plot(tiempo, Aixy, 'r-','linewidth', 1.5);
            %hold on
            %plot(tiempo, FS(:,6), 'b--','linewidth', 1.5);
            %xlim([-0.5 0.5])    

        else        % No Weygand
            
            % version all togeteher 2222
            %plot(tiempo(1:ps:end), FS(1:ps:end,6)-IS(1:ps:end,3), simbolo{3+i},'linewidth', 1.5);
            plot(tiempo, FS(:,6)-IS(:,3), simbolo{3+i},'linewidth', 1);
            if i == length(list) && tipo == 2
                legend('H&L', '30 W', '30', '60 W' ,'60', '120 W', '120' ,'location', 'southwest')
                legend box off
            end
            


            
            % version of the three graph separated  1111
            %plot(tiempo, FS(:,6)-IS(:,3), 'k-','linewidth', 1.5);
            %legend('H&L','Weygand', 'No Weygand', 'location', 'southwest')   
            %legend box off
            %set(gca,'FontWeight','bold');
            %set(gca,'linewidth',1)
            %xlabel('X (\mum)'); ylabel('\sigma (MPa)'); 
            %hold off
            %grid on

        end
        
        

    end
    
        %savefigure(list(i).name, 2)

end

disp FINISHED








%%%%%%%%%%%%  FUNCTIONS  %%%%%%%%%%%%%%%


function MakeMaps3(SigmaXX,SigmaXY, y,G,h)
figure(h)
 
    if G==1
    %plot(y,SigmaXX,'k-')
    %hold on
    plot(y,SigmaXY,'r-')%, 'linewidth', 1.5)
    %legend('W \sigma_{XX}','W \sigma_{XY}')
    end
    if G==2
    hold on
    plot(y,SigmaXX,'k--')
    plot(y,SigmaXY,'b--')
    %legend('\sigma_{XX}','\sigma_{XY}')
    end
    
    %hold off
    legend box off
    grid on
    set(gca,'FontWeight','bold');
    set(gca,'linewidth',1)
    xlabel('X (\mum)'); ylabel('\sigma (MPa)'); 
    xlim([-0.5 0.5]);
end

%%%-----------------------------------------------------------------------
% Adapt the 1d arrar of the stress components coming from Numodsi coupling
% to a six component column array

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
    plot(y,SigmaXX,'k-')%,'linewidth', 2)
    hold on
    plot(y,SigmaYY,'b-')%,'linewidth', 2)
    plot(y,SigmaXY,'r-')%,'linewidth', 2)
    plot(y,SigmaZZ,'g-')%,'linewidth', 2)
    legend('\sigma_{XX}', '\sigma_{YY}','\sigma_{XY}','\sigma_{ZZ}')
    end
    if G==2
    hold on
    plot(y,SigmaXX,'k.')
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
    plot(y,SigmaXX,'k--')
    plot(y,SigmaXY,'b--')
    %legend('\sigma_{XX}','\sigma_{XY}')
    end
    
    %hold off
    legend box off
    grid on
    set(gca,'FontWeight','bold');
    set(gca,'linewidth',1)
    xlabel('X (\mum)'); %ylabel('\sigma (MPa)'); 
    %xlim([-0.5 0.5]);
end
%%%%-----------------------------------------------------------------------

function out_Djon = Djon(Theo, Exp)
    if length(Theo) == length(Exp)
        out_Djon = (Exp-Theo)./Theo;
    else
        disp('Observed and predicted values do not have same size...')
        disp('Aborting D_Jon function... ')
    end
end%-----------------------------------------------------------------------


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

    nombre = este(1:length(este)-4);
    
    if n == 1
       fig_name = strcat('WAiry',num2str(nombre))

    else
       fig_name = strcat('MeshEffect',num2str(nombre));

    end

     print(fig_name,'-djpeg','-r1000'); 
end