clc
close all
clear all

%A = importdata('resFEMerateStrainControl/FEMLOAD.txt');
%B = importdata('resFEMerateStressControl/FEMLOAD.txt');
%D = importdata('resFEMerateStressControlrefined/FEMLOAD.txt');
%E = importdata('resNumodis/SIGEPS.csv');
F = importdata('../../res/FEMLOAD.txt');

%  "Time" "Real Time" 
%  "Sigma[11]($3)" "Sigma[22]($4)" "Sigma[33]($5)" "Sigma[12]($6)" "Sigma[13]($7)" "Sigma[23]($8)" 
%  "E[11]($9)" "E[22]($10)" "E[33]($11)" "E[12]($12)" "E[13]($13)" "E[23]($14)" 
%  "Ep[11]($15)" "Ep[22]($16)" "Ep[33]($17)" "Ep[12]($18)" "Ep[13]($19)" "Ep[23]($20)"


%
%   time[1]  Sigma33NC[2]   Sigma33C[3]  Str33NC[4]  Str33C[5]  Pstr[6]  RealStr[7]

%plot(A.data(:,7),D.data(:,3),'r*-'); xlabel('Strain'); ylabel('Stress (MPa)'), grid on
h1 = figure;
set(h1,'Position',[100 100 800 800])
subplot(2,2,1)
hold on
%plot(F.data(:,4),F.data(:,2),'r-')
%plot(F.data(:,5),F.data(:,3),'g-')
plot(F.data(:,7),F.data(:,3),'ko-'),xlabel('Strain'); ylabel('Stress (MPa)'), grid on
plot(F.data(:,4),F.data(:,3),'bo'),xlabel('Strain'); ylabel('Stress (MPa)'), grid on

subplot(2,2,2)
%figure
plot(F.data(:,1),F.data(:,7),'ro-'),xlabel('time'); ylabel('Strain'), grid on

subplot(2,2,3)
%figure
plot(F.data(:,1),F.data(:,2)-F.data(:,3),'ro-'),xlabel('time'); ylabel('Stress diff (Nc-C)'), grid on

subplot(2,2,4)
plot(F.data(:,1),F.data(:,6),'ro-'),xlabel('time'); ylabel('Ep'), grid on


%figure
%plot(E.data(:,11),E.data(:,5), 'G-')


