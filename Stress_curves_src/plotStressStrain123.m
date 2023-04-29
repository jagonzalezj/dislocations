clc
close all
clear all


format long
% % Material properties: Copper (Numodis)
E = 111216;            % Young Modulus  [Mpa]
% psi = 0.324;             % Poisson ratio
% b = 0.2552;              % Burgers vector [nm]
% miu = 42000;             % [MPa]
% %E = 2*miu*(1+psi);
% b = b * 0.001;     % nm    ==>  micrometer

InitStress = -0.0;    % MPa
     erate = -1e-9;   % /ns
       Ndt = 0.5;   % 0.05 ns
         N = 7;      % number of timesteps
        Ep = 7.3934672448298876E-009;      % plastic strain

         
   SzzNC = E*(erate*Ndt*N) + InitStress;
    SzzC = E*(erate*Ndt*N - Ep) + InitStress;
StrainNC = SzzNC/E;
 StrainC = SzzC/E;

disp(['Sigma  NC = ',num2str(SzzNC),''])
disp(['Sigma   C = ',num2str(SzzC),''])
disp(['StrainNCC = ',num2str(StrainNC),''])
disp(['Strain  C = ',num2str(StrainC),''])



% E = importdata('res/SIGEPS.csv');
% 
% % %  "Time" "Real Time" 
% % %  "Sigma[11]($3)" "Sigma[22]($4)" "Sigma[33]($5)" "Sigma[12]($6)" "Sigma[13]($7)" "Sigma[23]($8)" 
% % %  "E[11]($9)" "E[22]($10)" "E[33]($11)" "E[12]($12)" "E[13]($13)" "E[23]($14)" 
% % %  "Ep[11]($15)" "Ep[22]($16)" "Ep[33]($17)" "Ep[12]($18)" "Ep[13]($19)" "Ep[23]($20)"
% 
% h1 = figure;
% set(h1,'Position',[100 100 800 800])
% suptitle('SIGEPS NEW NUMODIS')
% subplot(2,2,1)
% hold on
% plot(E.data(:,2),E.data(:,17)*(-1),'k-', 'linewidth', 3),xlabel('time(ns)'); ylabel('plastic strain 33'), grid on
% plot(E.data(:,2),E.data(:,11)*(-1),'r-', 'linewidth', 3),xlabel('time(ns)'); ylabel('strain33'), grid on
% legend('Ep[33]','E[33]')
% subplot(2,2,2)
% plot(E.data(:,2),E.data(:,5)*(-1),'g-', 'linewidth', 3),xlabel('time (ns)'); ylabel('Stress33 (Mpa)'), grid on
% subplot(2,2,3)
% plot(E.data(:,11)*(-1),E.data(:,5)*(-1),'b-', 'linewidth', 3),xlabel('strain33'); ylabel('Stress33 (Mpa)'), grid on
%  subplot(2,2,4)
% plot(E.data(:,17)*(-1),E.data(:,5)*(-1),'b-', 'linewidth', 3),xlabel('Plastic strain33'); ylabel('Stress33 (Mpa)'), grid on
% 


F = importdata('res/FEMLOAD.txt');
% % % 
% % % %   step[1] time[2]  Sigma33NC[3]   Sigma33C[4]  Str33NC[5]  Str33C[6]
% % % %   Pstr[7]  ElmerStrain[8] ElmerStress[9]
% % % %
h1 = figure;
set(h1,'Position',[100 100 800 800])
suptitle('StressControl -1e-8(/ns)')

subplot(2,2,1)
plot(F.data(:,2),-F.data(:,5),'r-');xlabel('time (ns)'); ylabel('Strain')
hold on
plot(F.data(:,2),-F.data(:,6),'k-')
plot(F.data(:,2),F.data(:,8),'b-');xlabel('time (ns)'); ylabel('Strain Elmer')
legend('NC Strain','C Strain', 'ElmerStrain', 'location', 'best')
grid on

subplot(2,2,2)
plot(F.data(:,2), -F.data(:,7),'k-');xlabel('time (ns)'); ylabel('plastic strain')
grid on

subplot(2,2,3)
plot(F.data(:,2),-F.data(:,3),'r-');xlabel('time (ns)'); ylabel('Stress')
hold on
plot(F.data(:,2),-F.data(:,4),'k-')
%plot(F.data(:,1),F.data(:,7)*E,'b-')
plot(F.data(:,2),F.data(:,9),'b-')
legend('NC Stress','C Stress','ElmerStress', 'location', 'best')
grid on

subplot(2,2,4)
plot(-F.data(:,5),-F.data(:,3),'r-');xlabel('Strain'); ylabel('Stress')
hold on
plot(-F.data(:,5),-F.data(:,4),'k-')
plot(-F.data(:,5),F.data(:,9),'b-')
legend('NC Stress','C Stress','ElmerStress', 'location', 'best')
grid on

