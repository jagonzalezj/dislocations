%%%------------------------------------------------------------------------
% This code was created to read an spreadSheet from paraview with the
% outputs of ElmerNumodis coupling and adapt the point coordinates to the
% ones that are read by to coupling to export the self stress.
% It makes the positions as: 
%
%  !   total# x1 y1 z1 x2 y2 z2 ....
%
% also the stress coming from the same spreadsheet will be acomodate using
% the format of output of numodis (the spreadsheet should have all the data
%   points, displacement, strain and stress)
%
%  ! Spreasheet  Elmer : XX XY XZ YY YZ ZZ
%             ! Numodis: XX YY ZZ XY XZ YZ
%
%  outputs:
%              pos.txt ===> Positions as Numodis
%           Stress.txt ===>  Stresses as Numodis
 %%%------------------------------------------------------------------------

clc
close all 
clear all

% Spreadsheet from Paraview
A = importdata('../SelfStressInf2/FEMBoundLineRefWeygand.csv'); 
  
% Makes Pos and Stress one array vector
[Pos, Stress] = EnlargePos(A);


disp("Finished ExtractAdapt")



function [Pos, Stress] = EnlargePos(A)
       Pos = [];
    Stress = [];
    for i= 1:length(A.data(:,1))

        Pos(1,3*i-2) = A.data(i,1);  % Positions
        Pos(1,3*i-1) = A.data(i,2);  
        Pos(1,3*i)   = A.data(i,3);  

        Stress(1,6*i-5) = A.data(i,15);  % stresses
        Stress(1,6*i-4) = A.data(i,18);  
        Stress(1,6*i-3) = A.data(i,20);  
        Stress(1,6*i-2) = A.data(i,16);  
        Stress(1,6*i-1) = A.data(i,17);  
        Stress(1,6*i)   = A.data(i,19); 
    end

    fid1=fopen('pos.txt','w');
    fprintf(fid1, '%d \n', length(Pos));
    fprintf(fid1, '%f ', Pos');
    fclose(fid1);


    fid2=fopen('StressParaview.txt','w');
    fprintf(fid2, '%f ', Stress');
    fclose(fid2);
end




