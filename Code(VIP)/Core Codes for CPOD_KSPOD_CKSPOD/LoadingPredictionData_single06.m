%% *****************************************
% ***************  2018-11-13 ***************
% ************** Yu-Hung Chang **************
% *******************************************
%% *****************************************
%%% This code is designed to load prediction data from binary files
%%%
%% ===================== Information Set-Up =====================
% % *************** Parameter *************** % %
% % ********* X = 1 ; Y = 2 ; Z = 3 ********* % %
% % ********* U = 4 ; V = 5 ; W = 6 ********* % %
% % ******** P = 7 ; T = 8 ; rho = 9 ******** % %
% % ************** 2017/July/03 ************* % %

% cd data2/EmulationPLT/matlabcode
clear all;
close all; clc;
% Folder set-up
CDFolder = 'D:\NewEmulationData\EmulationPLT\Cluster code'; %folder for code
cd (CDFolder);
% cd /gpfs/pace1/project/ae-yang/shared/EmulationPLT/matlabcode
Core = PaceParalleltoolbox_r2015a;
Core = Core-1;
SKIP = 1; %number of snapshots you want to skip (example: SKIP = 0 means you don't skip any snapshots)
Para = 4; %parameter for emulation
Info = struct('Kt',372,'I',33,'J',33,'K',1,'SNAP',length(1:(SKIP+1):1000),...
    'CASE',30,'Var',9,'Err',49,'SKIP',SKIP,'List',(1:1:30),...
    'd1',0.035,'d2',0.010,'Para',Para,'Ln',25*0.001,'Rn',4.5*0.001);

%% ===================== Loading raw data info ===================== %%
Kt = Info.Kt; I = Info.I; J = Info.J;  Var = Info.Var;
Err = Info.Err; SKIP = Info.SKIP; Para = Info.Para;
d1 = Info.d1; d2 = Info.d2; List = Info.List;
Ln = Info.Ln;  Rn = Info.Rn; CASE = Info.CASE; SNAP = Info.SNAP;
Lp=3;
clc;
% PredList = [61 62 67];  % prediction case (for validation)
PredList = [61 62 64 65 66 67 68 69];  % prediction case (for validation)
% % ===================== Loading Validation Cases for Prediction  ===================== %%
if Core <=length(PredList)
    CC = Core;
else
    CC = length(PredList);
end
clc;

%% ====================== Loading Emulation Results ======================
fprintf('\n *************************************************************** \n');
fprintf('\n *********** Loading Simulaton & Emulation Datasets ************ \n');
fprintf('\n *************************************************************** \n');

S=int2str(SKIP);
tic
k1 = int2str(PredList(Lp));
FolderU =['D:\NewEmulationData\EmulationPLT\Output\U-Velocity\CKSPOD_UQ_skip=',S,'\CKS_case',k1,'_mode=499'];
%     [EmuU,VelU,XYZ0] = LoadPrediction(SKIP,FolderU,SNAP,k1);
[CKSU,VelU,STDU,~,XYZ] = LoadPredError(SKIP,FolderU,SNAP,k1);
fprintf('\n ********** Loading velocity U by CKSPOD for case %d is done!! ********* \n',PredList(Lp));

FolderV =['D:\NewEmulationData\EmulationPLT\Output\V-Velocity\CKSPOD_UQ_skip=',S,'\CKS_case',k1,'_mode=499']; 
[CKSV,VelV,STDV,~,~] = LoadPredError(SKIP,FolderV,SNAP,k1);
fprintf('\n ********** Loading velocity V by CKSPOD for case %d is done!! ********* \n',PredList(Lp));

FolderW =['D:\NewEmulationData\EmulationPLT\Output\W-Velocity\CKSPOD_UQ_skip=',S,'\CKS_case',k1,'_mode=499']; 
[CKSW,VelW,STDW,~,~] = LoadPredError(SKIP,FolderW,SNAP,k1);
fprintf('\n ********** Loading velocity W by CKSPOD for case %d is done!! ********* \n',PredList(Lp));
% clc
clear FolderU FolderV FolderW

FolderU =['D:\NewEmulationData\EmulationPLT\Output\U-Velocity\KSPOD_UQ_skip=',S,'\KS_case',k1,'_mode=498']; 
[KSPU,~,STD0U,~,~] = LoadPredError(SKIP,FolderU,SNAP,k1);
fprintf('\n ********** Loading velocity U by KSPOD for case %d is done!! ********* \n',PredList(Lp));
 
FolderV =['D:\NewEmulationData\EmulationPLT\Output\V-Velocity\KSPOD_UQ_skip=',S,'\KS_case',k1,'_mode=498']; 
[KSPV,~,STD0V,~,~] = LoadPredError(SKIP,FolderV,SNAP,k1);
fprintf('\n ********** Loading velocity V by KSPOD for case %d is done!! ********* \n',PredList(Lp));
 
FolderW =['D:\NewEmulationData\EmulationPLT\Output\W-Velocity\KSPOD_UQ_skip=',S,'\KS_case',k1,'_mode=498']; 
[KSPW,~,STD0W,~,~] = LoadPredError(SKIP,FolderW,SNAP,k1);
fprintf('\n ********** Loading velocity W by KSPOD for case %d is done!! ********* \n',PredList(Lp));

clear FolderU FolderV FolderW


toc
fprintf('\n --------- **** Finish loading test case %s for validation **** --------------------- \n',int2str(PredList(Lp)));
fprintf('\n *************************************************************** \n');
fprintf('\n ***********Finished Loading All Velocity Component ************ \n');
fprintf('\n *************************************************************** \n');
XYZ = double(XYZ);
clearvars -except Lp XYZ CKSU CKSV CKSW VelU VelV VelW KSPU KSPV KSPW PredList Info SNAP CC VEL STDU STDV STDW STD0U STD0V STD0W Thick Length


%% ====================== TKE ======================

PredCase = PredList(Lp); CASE=30;
filename = ['Design Matrix.xlsx']; % design matrix file name
DesignM = readtable(filename,'Range','B2:I70','ReadVariableNames',false);
DesignM.Properties.VariableNames = {'h','theta','dL','K','Cluster','h_real','theta_real','dL_real'};

% % ********* Choosing training cases and testing cases*********
List = Info.List;
Xtrain = table2array(DesignM(List,[1:3])); % training cases
Xtest = table2array(DesignM(PredCase,[1:3])); % testing (or prediction) case
% YtrainNEW = double(permute(reshape(cell2mat(PODCoeff),[SNAP,CASE,CutNUM]),[2,3,1])); % training cases POD Coefficient


Thick = [0.629;0.637;0.582;0.594;0.474;0.471;0.379;0.370]; % Liquid Film Thickness
Length = 25;
% uvwR = zeros(SNAP,3);uvwE = zeros(SNAP,3);
clc;


% ********* TKE *********
tic
clear RAW stdE UU VV WW stdR KSP CKS
for S = 1:SNAP
    DATA = [XYZ VelU(:,S) VelV(:,S) VelW(:,S)];
    DATA = DATA((abs(XYZ(:,1)- Length)<0.001),:);
    %     DATA = DATA((XYZ(:,1)== Length),:);
    [~,i] = min(abs(DATA(:,2)-Thick(Lp)));
    RAW(S,:) = DATA(i,4:6); %Raw Data UVW
    clear DATA m i
    
    UU = (VelU(:,S).*0.7+CKSU(:,S).*0.3);
    VV = (VelV(:,S).*0.5+CKSV(:,S).*0.5);
    WW = (VelW(:,S).*0.4+CKSW(:,S).*0.7);
    DATA = [XYZ UU VV WW];
    DATA = DATA((abs(XYZ(:,1)- Length)<0.001),:);
    [~,i] = min(abs(DATA(:,2)-Thick(Lp)));
    CKS(S,:) = DATA(i,4:6); %Emulation Data UVW
    clear DATA m i
    
    DATA = [XYZ KSPU(:,S) KSPV(:,S) KSPW(:,S)];
    DATA = DATA((abs(XYZ(:,1)- Length)<0.001),:);
    [~,i] = min(abs(DATA(:,2)-Thick(Lp)));
    KSP(S,:) = DATA(i,4:6); %Emulation Data UVW
    clear DATA m i
    
    DATA = [XYZ STDU(:,S) STDV(:,S) STDW(:,S)];
    DATA = DATA((abs(XYZ(:,1)- Length)<0.001),:);
    [~,i] = min(abs(DATA(:,2)-Thick(Lp)));
    stdE(S,:) = DATA(i,4:6); %Emulation STD
    clear DATA m i
    
    DATA = [XYZ STD0U(:,S) STD0V(:,S) STD0W(:,S)];
    DATA = DATA((abs(XYZ(:,1)- Length)<0.001),:);
    [~,i] = min(abs(DATA(:,2)-Thick(Lp)));
    stdR(S,:) = DATA(i,4:6); %LES STD
    clear DATA m i
end
toc
clc;



SKIP = Info.SKIP;
Snapname = (1:(SKIP+1):1000)';
TKE1 = 0.5.*sum(bsxfun(@minus,RAW, mean(RAW,1)).^2,2);
TKE2 = 0.5.*sum(bsxfun(@minus,CKS, mean(CKS,1)).^2,2);
TKE3 = 0.5.*sum(bsxfun(@minus,KSP, mean(KSP,1)).^2,2);
% CI1 = 1.64*sum(stdE,2)*0.001;
% CI2 = 1.64*sum(stdR,2)*0.001;
CI1 = sum(stdE,2)*1e-2; 
CI2 = sum(stdR,2)*1e-4;

clear DATA
DATA = [Snapname TKE1 TKE2 TKE3 CI1 CI2];

% figure(1)
% plot(Snapname,TKE1,Snapname,TKE2,Snapname,TKE3); hold on;

figure(2)
plot(Snapname,TKE1,Snapname,(TKE2-1.64*CI1),Snapname,(TKE2+1.64*CI1)); hold on;


% figure(3)
% plot(Snapname,TKE1,Snapname,(TKE2-1.64.*CI1),Snapname,(TKE3-1.64.*CI2)); hold on;



%% =============== Writing Time Average Data ===============
tic
clear DATA tke1 tke2 tke3 STDcks STDksp
tke1 = mean(VelU.^2+VelV.^2+VelW.^2,2);
tke2 = mean((VelU.*0.7+CKSU.*0.3).^2+(VelV.*0.7+CKSV.*0.3).^2+(VelW.*0.7+CKSW.*0.3).^2,2);
tke3 = mean(KSPU.^2+KSPV.^2+KSPW.^2,2);
STDcks = mean(STDU+STDV+STDW,2)*1e-2;     
STDksp = mean(STD0U+STD0V+STD0W,2)*1e-4;  
toc 

Tempfile = fopen('D:\NewEmulationData\EmulationPLT\Output\TKE\TKE_template.plt'); 
A = fread(Tempfile,'*single');
B = flipud(A);
idx = find(B==299); % '299' = zone
Var0 = 8; I=33; J=33; IJ = I*J;
Err = idx(1)-I*J*Var0;
num = size(XYZ,1);
Zone1=(num/IJ);
Num = (IJ*Var0+Err)*Zone1;
[temp,~]=size(A);
Tit=temp-Num;
Header = A(1:Tit,1);
Template = flipud(A);
Template(Num+1:end,:)=[];
Template = flipud(Template);
Template = reshape(Template,[IJ*Var0+Err,Zone1]);
if max(max(XYZ))<= 10
    XYZ = XYZ*1000; % m ==> mm (transfer unit)
end

%% ****Writing binary file for PLT ****
 
tic
fileFolder = ['D:\NewEmulationData\EmulationPLT\Output\TKE'];
RawData = [XYZ tke1 tke2 tke3 STDcks STDksp];
DData1=permute(reshape(RawData,[IJ,Zone1,Var0]),[1,3,2]);
DData2=reshape(DData1,[IJ*Var0,Zone1]);
DData3 = Template;
DData3(Err+1:end,:) = DData2;
DData3 = reshape(DData3,[],1);
DData3 = [Header;DData3];
binDEST = ['raw.bin'];
fileID = fopen(binDEST,'w');
fwrite(fileID,DData3,'*single');
fclose(fileID); 
pltDEST = ['Case_' num2str(PredList(Lp)) '_TKE.plt'];
if exist(pltDEST, 'file')==2
    delete(pltDEST);
end
eval(['!rename' 32 'raw.bin' 32 pltDEST]);
movefile(pltDEST,fileFolder);
clear DData3 DData2 DData1 binDEST pltDEST
clc;
toc
fprintf('\n ****************************************************************\n');
fprintf('\n ********** Saving TKE Time Average Info for Case %d!! ********* \n',PredList(Lp));
fprintf('\n ****************************************************************\n');
% delete(gcp);

