%%%------------------------------------------------------------------------
%   Plot the Difference between the "total nodal stress" at a surface of the 
%   domain used in the coupling and the "Self Stress" at the same position.
%   
%   After a coupling timestep the nodes positions are moved by the stress
%   of the dislocations. It means the the matlab fucntions "repeat" and
%   "unic" will not work properly.
%   To correct it a "Ficticious" node postion array is created to make
%   psosible the use of "Contourn" mapping fucntions and "Meshgrids"
%
%   Working flow:
%    1. Run the Coupling and extract a planar slice saved as FEMstressSurface.txt
%    2. Run ExtractAdapt.m and create the files pos.txt (position as
%       Numodis) and Stress.txt (Stress as numodis)
%    3. Run the coupling one more time now readin the new pos.txt file and
%       creating the new file Selfstress.txt and/or WeygandStress.txt
%    4. Finally run this code and compare the SelfStress.txt to the
%       Stress.txt. Then to plot the diference a trick should be made in
%       order of make an efficient meshgrid. 
%
%%%------------------------------------------------------------------------

clc
close all 
clear all

Coord = importdata('pos.txt'); % import one file data (numodis format)
Coord = Coord(2:end);         % avoid first number (total points)

% create the 3d vectors to plot.
for i = 1: length(Coord)/3
    X(i)=[Coord(3*i - 2)]/10000;
    Y(i)=[Coord(3*i - 1)]/10000;
    Z(i)=[Coord(3*i)]/10000;
end

% Optional : plot the plane where stress is requested
%plot3(X,Y,Z,"k."); xlabel("X(1)"); ylabel("Y(2)"); zlabel("Z(3)"); grid on; axis square

%Here prepare the plotting functions
repX = unique(X);
repY = unique(Y);
repZ = unique(Z);

%SelfStress = load('SelfStress.txt');
StressPar = load('StressParaview.txt');

if length(repX)== 1
    disp('ploting normal to YZ')
   
elseif length(repY) == 1
    disp('ploting normal to XZ')
    [arr refarr] = updateposs(Coord, 30, 30, [-4000.0 4000.0], [-75000 75000], 2, repY);
    arr
    %sortrows(arr,[3,1])
    %xp= reshape(arr(:,1), [30,30])
    
elseif length(repZ)==1
    disp('ploting normal to XY')
end


% Separate the six components of the SelfStress make it a fucntion
% m = 0;
%     for i = 1 : length(al)
%         for j = 1 : length(be)
%             m= m + 1;
%             Sxx(i,j) = Stress(6*m-5); 
%             Syy(i,j) = Stress(6*m-4);
%             Szz(i,j) = Stress(6*m-3);
%             Sxy(i,j) = Stress(6*m-2);
%             Sxz(i,j) = Stress(6*m-1);
%             Syz(i,j) = Stress(6*m);
%         end
%     end 
    
    
    
    
    %%%%%%%%%     FUNCTIONS     %%%%%%%%%%%
   
    function [NewArr,ArtArr] = updateposs(arr,step1, step2, lim1, lim2, eje, rep);
    % arr is the pos array from spreadsheet of paraview
    % step1 and 2 refer to step of mesh generated -40:step:40
    % lim1= [-40,40] ....  same for lim2
    % 1 and 2 refer to X and Y respectivelly.
     
        pass1 = (abs(lim1(1))+abs(lim1(2)))/(step1-1);  % x step
        pass2 = (abs(lim2(1))+abs(lim2(2)))/(step2-1);  % y step

        a = lim1(1): pass1: lim1(2);   % X vector
        b = lim2(1): pass2: lim2(2);   % Y vector
    
        % create new array to pass the new poss values
        % convert [1d] pos array to [n,3] x,y,z array
        NewArr(length(arr)/3, 3) = [0]; 
        for j = 1:length(arr)/3
            NewArr(j,:)=[arr(3*j-2) arr(3*j-1) arr(3*j)] ;   
        end       
    
        % creating the new artifical array [n,3] x,y,z to be compared
        if length(arr)/3 == length(a)*length(b)
            ArtArr(length(a)*length(b),3) = [0]; % artifical array
            xflag = 1;
            yflag = 1;
            zflag = 1;         
                switch eje  % select "Axis normal" case      
                    
                    case 2  % case when AXIS is <Y> normal
                        for k = 1:length(a)*length(b) % loop to fill  ArtArr
                            ArtArr(k,:)=[a(xflag) rep b(zflag)];
                            if zflag == length(b) && xflag < length(a) 
                                xflag = xflag + 1;
                                zflag = 0;
                            end
                            zflag = zflag + 1;
                        end      
                    
                    case 3  % case when AXIS is <Z> normal
                        for k = 1:length(a)*length(b) % loop to fill  ArtAr
                            ArtArr(k,:)=[a(xflag) b(yflag) rep];
                            if yflag == length(b) && xflag < length(a) 
                                xflag = xflag + 1;
                                yflag = 0;
                            end
                            yflag = yflag + 1;
                        end
                end      
        else
            disp("the length of both arrays are not the same")
        end
     end