%% *****************************************
% ***************  2018-11-13 ***************
% ************** Yu-Hung Chang **************
% *******************************************
%% *****************************************

close all;clear all; clc; % close all & clean all & clear all in the command window
CDFolder = 'D:\NewEmulationData\EmulationPLT\Cluster code'; %name of the directory for code 
cd (CDFolder); % set the directory

Main_Setting; % code for Main setting


%% ******** Code for loading training data reduction (choose CPOD, KSPOD, or CKSPOD) *******
clc;
Main_TrainingData; % code for loading training data 

 

%% ******** Code for data reduction (choose CPOD, KSPOD, or CKSPOD) *******
clc;
% Main_DataReduction_CKSPOD;
Main_DataReduction_KSPOD;
% Main_DataReduction_CPOD;
 
%% ******** Code for GP fit and validation (choose CPOD, KSPOD, or CKSPOD) *******
clc; 
% Main_GPfit_validation_CKSPOD; 
% Main_GPfit_validation_KSPOD; 
% Main_WritingPODInfo;




