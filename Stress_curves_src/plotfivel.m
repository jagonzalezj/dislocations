clc 
close all
clear all

%A = load('FigureFivel.tsv');
A = load('newloopfivel.tsv');


top = max(A(:,3))

plot(top-A(:,3),A(:,6)*1e-6,'r-')
hold on
plot(top-A(:,3),A(:,7)*1e-6,'b-')
plot(top-A(:,3),A(:,8)*1e-6,'g-')
plot(top-A(:,3),A(:,9)*1e-6,'y-')
plot(top-A(:,3),A(:,10)*1e-6,'k-')
plot(top-A(:,3),A(:,11)*1e-6,'m-')

ylabel("Stress (MPa)")
xlabel("Distance from the top (A)")
xlim([0.0 20000])
legend('Sxx',	'Sxy',	'Sxz',	'Syy',	'Syz', 'Szz')