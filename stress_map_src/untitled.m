clc
clear all
close all


% Material properties: Aluminum
% E = 70000 ;            % Young Modulus  [Mpa]
% psi = 0.35;            % Poisson ratio
% b = 0.285;             % Burgers vector [nm]
% miu = E/(2*(1+psi));   % Shear Modulus  [Units of E]

% % Material properties: Copper (Numodis)
%E = 111216 ;            % Young Modulus  [Mpa]
psi = 0.324;             % Poisson ratio
b = 0.2552;              % Burgers vector [nm]
%miu = E/(2*(1+psi));    % Shear Modulus  [Units of E]
miu = 42000;             % [MPa]
E = 2*miu*(1+psi);

% % % Material properties: Copper (Liu 2009)
% %E =  ;            % Young Modulus  [Mpa]
% psi = 0.347;            % Poisson ratio
% b = 0.256;             % Burgers vector [nm]
% %miu = E/(2*(1+psi));   % Shear Modulus  [Units of E]
% miu = 42000;

% units conversions if needed
b = b * 0.001;     % nm    ==>  micrometer
prim = (miu*b)/(2*pi*(1-psi));

%%%%----------------------------------------------------------------------
%%%%%%   Calculation of image stress (Theory from Hirt and Lothe)
%%%%%%   There is:
%%%%%%       Total_Stress 
%%%%%%       (Total_Stress) - Self_Stress = Image_Stress
%%%%-----------------------------------------------------------------------
%%% Acomodate the distance to surface to mathc with self stress plot
%%% the "0" is at the boundarie, so change it :)
X = [-8:0.05:0];  % micrometers
Y = [-7:0.05:7];
[x,y]=meshgrid(X,Y);
%%%%-----------------------------------------------------------------------
%%%%%%      Total Stress Edge (Hirt and Lothe {Singular}) 
%%%%-----------------------------------------------------------------------
l = 1;
r = sqrt((x-l).^2 + y.^2);
I_SigmaXX = -prim * (4*l*x.*y)./(r.^6) .* (3.*(l-x).^2 - y.^2);
I_SigmaYY =  prim * (2*l)./(r.^6) .* ( 4*y.*(l-x).^3 + 6*x.*y.*(l-x).^2 + 4*(y.^3).*(l-x) - 2*x.*y.^3);
I_SigmaXY = -prim * (2*l)./(r.^6) .* ( (l-x).^4 + 2*x.*(l-x).^3 - 6*x.*(y.^2).*(l-x) - y.^4);
I_SigmaZZ =  prim * (8*l*psi)./(r.^6) .* (y.*(l-x).^3 + (l-x).*(y).^3); 
%%%%-----------------------------------------------------------------------

%%%%-----------------------------------------------------------------------

%--% PLOTIING Total_Stress 
%%%%-----------------------------------------------------------------------
h2=figure(2);
subplot(2,2,2); contourf(x, y, I_SigmaXX,50, 'ShowText', 'off', 'LineColor', 'none'); %hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('I_SigmaXX (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal; 
colorbar
%caxis([-3.8 3.8])

subplot(2,2,3); contourf(x, y, I_SigmaYY,20, 'ShowText', 'off', 'LineColor', 'none'); %hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('I_SigmaYY (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal; 
colorbar

subplot(2,2,1); contourf(x, y, I_SigmaZZ, 20, 'ShowText', 'off', 'LineColor', 'none'); %hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('I_SigmaZZ (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mum)'); 
axis equal; 
colorbar
%caxis([-3.3 3.3])

subplot(2,2,4); contourf(x, y, I_SigmaXY,20, 'ShowText', 'off', 'LineColor', 'none'); %hold on; plot(l,0,'k.','MarkerSize',30, 'linewidth', 10)
title('I_SigmaXY (Mpa): Inf. Edge dislocation'); xlabel('X (\mum)'); ylabel('y (\mumm)'); 
axis equal; 
colorbar
%caxis([-0.01 0.01])
colormap(jet)
set(h2,'position',[10 10 1000 1000]); 


