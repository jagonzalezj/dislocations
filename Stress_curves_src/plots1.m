clc
close all
clear all

format short 
% Nucleation Constants definititions
Ncites=1;
Tm = 700.0;              % K
nu_0 = 3.14e11;          % s^-1
Kb = 1.380649e-23 ;      % J/K --> (m^2 kg)/(s^2 K) 
T =650;                  % K
numodisdt = 0.05;      % (ns)
alpha = 1.0;

D_E = [20 0.2 20] * 1.60218e-19  % convert to joules

% set of equations    
N= length(D_E);
Nu = Ncites*nu_0 * exp (-(D_E *(1-T/Tm))/(Kb*T));
dt = 1/(alpha*sum(Nu));
D_N = Nu * dt;
D_N( length(D_N)+1 ) = 1 - sum(D_N); % include the "Do nothing"

% calculate nucleation timesteps
NSTEP = (numodisdt/(dt*1E9));

% generate cinturon
Cint = [0.0 cumsum(D_N)];

% uniform randon number generation
randomnum = unifrnd(0,1,1000,1);
%figure  %%%%%%% fist graph
% plotforpaper
% hist(randomnum,10)
% set(gca,'FontWeight','bold');
% set(gca,'linewidth',1);
% xlabel('Probability'); ylabel('Counts');
% stop

% statistics bar for jonathan
B(length(Cint)-1) =  zeros;


disp(['Numodis dt: [',num2str(numodisdt),']'])
disp(['DE: [',num2str(D_E/1.60218e-19),']'])
disp(['nKMC: ',num2str(NSTEP),''])
disp(['1/nKMC : ',num2str(1/NSTEP),''])
disp(['Calculated dt: ',num2str(dt),' (s)   ',num2str(dt*1E9), '(ns)' ])
disp(['DN: [',num2str(D_N),']'])
disp(['Markov: [',num2str(Cint),']'])

plotforpaper
yyaxis right
%hAx=gca;
% here is important
for i =1 : length(randomnum) ;% loop over the randomnumber    
    for j = 1 : length(Cint)-1; % loop over steps in Cint
        if j == length(Cint)-1;
            %disp("Do nothing case")
            plot(randomnum(i),rand(1),'k.','markersize',8)
             B(j) = B(j)+1;
           break 
        end
        
        if randomnum(i)<Cint(j+1) && randomnum(i)>Cint(j);
           plot(randomnum(i),rand(1),'k.','markersize',8)
           hold on
           %fprintf('nucleate case:%i %f entre  %f y %f \n',j, D_N(j), Cint(j), Cint(j+1));
           B(j) = B(j)+1;
           break
        end        
    end
end
xlim([0 1]);
%hAx.LineWidth=5;

% another Figure
for  h = 1:length(Cint);
   xline(Cint(h),'k-', 'linewidth', 2);
end

%Draw rectangles
yyaxis left
hold on
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4660 0.6740 0.1880]};
for k =1:length(Cint)-1;
 rectangle('Position', [Cint(k), 0, Cint(k+1)-Cint(k), B(k)*1], 'EdgeColor', 'none', 'FaceColor', colors{k}, 'LineWidth', 4);
end

set(gca,'FontWeight','bold');
set(gca,'linewidth',1);
xlabel('Cumulative Probability')
ylabel('Counts')

%set(gca,'YAxisLocation','right','ytick',[]);


% set(gca,'ytick',[]);
% yLimits = get(gca,'YLim');
% rectangle('Position', [0, 0, 1, yLimits(2)], 'LineWidth', 3);
% pbaspect([10 1 1]);
% yLimits = get(gca,'YLim');
% 
% 

% %%%% PAPER FIGURE
% figure
% plotforpaper;
% X = categorical({'250','650'});
% y = [ 0.49997  0.37284;   0.49997 0.37284; 0.0 0.25432];
% h=bar(X,y','stacked','EdgeColor','k', 'linewidth', 1.2);
% set(h, {'DisplayName'}, {'c1','c2','s'}');
% legend();
% set(gca,'FontWeight','bold');
% set(gca,'linewidth',1);
% xlabel('Temperature (K)'); ylabel('Probability');




function plotforpaper();
    fig = figure;
    fig.Units  = 'centimeters';
    fig.Position(3) = 9;
    fig.Position(4) = 9;
    set(fig.Children,'FontName','Times','FontSize',10);
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));% remove white space
    fig.PaperPositionMode   = 'auto';
    axis square;
    % %set(gcf, 'color', 'none'); set(gca, 'color', 'none');
end %----------





% hAx=gca;                    % create an axes
% hL=plot(randn(5,1));        % and a line underneath
% hAx.LineWidth=2;            % set the axis linewidth for box/ticks
% hLg=legend('A');            % put a legend on the figure
% hLg.LineWidth=1;            % make the legend axes box linewidth smaller