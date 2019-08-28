%% *****************************************
% ***************  2018-11-13 ***************
% ************** Yu-Hung Chang **************
% *******************************************
%% *****************************************
%% ========================== Information Set-Up ===========================
% % ***************************** Parameter ************************************* % %
% % ********* # represent the order of each parameter in a raw data file ******** % %
% % *************************** X = 1 ; Y = 2 ; Z = 3 *************************** % %
% % *************************** U = 4 ; V = 5 ; W = 6 *************************** % %
% % ************************** P = 7 ; T = 8 ; rho = 9 ************************** % %
% % ***************************************************************************** % %

% cd data2/EmulationPLT/matlabcode
% clear all; 
close all; clc;
% Folder set-up
CDFolder = 'D:\NewEmulationData\EmulationPLT\Cluster code'; %folder for code 
cd (CDFolder); 
% cd /gpfs/pace1/project/ae-yang/shared/EmulationPLT/matlabcode
Core = PaceParalleltoolbox_r2015a;
Core = Core-1;
SKIP = 2; %number of snapshots you want to skip (example: SKIP = 0 means you don't skip any snapshots)
Para = 7; %parameter for emulation
Info = struct('Kt',372,'I',33,'J',33,'K',1,'SNAP',length(1:(SKIP+1):1000),...
    'CASE',30,'Var',9,'Err',49,'SKIP',SKIP,'List',(1:1:30),...
    'd1',0.035,'d2',0.010,'Para',Para,'Ln',25*0.001,'Rn',4.5*0.001);
% ======================================= About the "Info" =======================================
% % ********* KT = number of zones; I & J are the I and J For Tecplot binary files ************ % %
% % ********* CASE = total numbers of training cases ****************************************** % %
% % ********* SNAP = total numbers of snapshots *********************************************** % %
% % ********* Var = total variables in Tecplot binary files *********************************** % %
% % ********* Err = total lines of errors reading by MATLAB from Tecplot binary files ********* % %
% % ********* SKIP = the number of snapshots you want to skip during reading ****************** % %
% % ********* Ln = leangth of the swirl injector (m) ****************************************** % %
% % ********* Rn = radius of the swirl injector (m) ******************************************* % %
% % ********* d1 = the cut-off range along with X (m) ***************************************** % %
% % ********* d2 = the cut-off range along with Y (m) ***************************************** % %
