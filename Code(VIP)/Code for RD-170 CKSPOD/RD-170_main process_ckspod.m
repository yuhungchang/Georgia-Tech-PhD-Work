
cd D:\NewEmulationData\RD-170Group\code
%% ===================== Loading raw data ===================== %%
clear all; close all; clc;
Core = PaceParalleltoolbox_r2015a;
% Ln = 25*0.001;  Rn = 4.5*0.001;
% % *************** Parameter *************** % %
% % ********* X = 1 ; Y = 2 ; Z = 3 ********* % %
% % ********* U = 4 ; V = 5 ; W = 6 ********* % %
% % ******** P = 7 ; T = 8 ; rho = 9 ******** % %
% % ***************************************** % %
Para = 10;
I=33;
J=17;
% Var = 23;
d1 = 0.175;
d2 = 7.8e-2;
% d2 = 8.45e-2;
% d2 = 8.75e-2;
d3 = 0.025;
r1 = 0.1035;
r2 = 5.62e-3;
r3 = 6.365e-3;

Range = [d1 d2 d3 r1 r2 r3];
SNAP = 450; % length(1:(SKIP+1):1000);
CASE = 10;
DataD = cell(CASE,1);
DataU = cell(CASE,1);
DataR = cell(CASE,1);
DataUB = cell(CASE,1);
xyzD = cell(CASE,1);
xyzU = cell(CASE,1);
xyzR = cell(CASE,1);
xyzUB = cell(CASE,1);
fprintf('\n *********** Loading raw data *********** \n');
parfor (Lp = 1:CASE)  %Case number
    k1 = int2str(Lp);
    DirOut = ['D:\NewEmulationData\RD-170Group\Case',k1,''];
    if (Lp ==3 | Lp==6)
        Var = 24;
    else
        Var = 23;
    end
    tic
    S = 1;
    for fil = 1:SNAP
        k2 = int2str(fil+10);
        filename = ['',k2,'.plt'];
        fileDEST  = fullfile(DirOut,filename);
        [xyzD{Lp},xyzU{Lp},xyzR{Lp},xyzUB{Lp},DataD{Lp}(:,S),DataU{Lp}(:,S),DataR{Lp}(:,S),DataUB{Lp}(:,S)] = PLTloaderRD170(I,J,Var,Para,Range,fileDEST);
        S= S+1;
    end
    T(Lp) = toc;
    DataD{Lp} = double(DataD{Lp});
    xyzD{Lp} = double(xyzD{Lp});
    
    DataU{Lp} = double(DataU{Lp});
    xyzU{Lp} = double(xyzU{Lp});
    
    DataR{Lp} = double(DataR{Lp});
    xyzR{Lp} = double(xyzR{Lp});
    clc
    fprintf('\n *********** Complete loading Case %d for %d files, Time = %5.2f sec *********** \n',Lp,S,T(Lp));
end

S = size(DataD{1},2);
fprintf('\n *********** Complete loading %d cases for %d snapshots *********** \n',CASE,S);

% 139-212//215-232//235-252//255-272//275-292
%% ===================== Common Grid  ===================== %%
clc;
delimiterIn = ' ';
headerlinesIn = 30;

FileDir= 'D:\NewEmulationData\RD-170Group';
MeshFile = {'Mesh_zone57-96_merge.dat';'Mesh_zone97-138_merge.dat';'Mesh_zone139-292_select.dat';'Mesh_zone139-294_merge.dat'};

% [Common,Label,commonXYZ,Node] = CommonRD170(CASE,SNAP,xyzU,xyzD,DataU,DataD,MeshFile,FileDir,delimiterIn,headerlinesIn);
% [Common,commonXYZ,Node] = CommonRD170down(CASE,SNAP,xyzD,DataD,MeshFile,FileDir,delimiterIn,headerlinesIn);
[Common,commonXYZ,Label,Node] = CommonRD170(CASE,SNAP,xyzU,xyzD,DataU,DataD,MeshFile,FileDir,delimiterIn,headerlinesIn);



%% ===================== Main Part (POD) ===================== %%
PredList = [7,10];
for PRED = 1:length(PredList)
    clearvars -except Common Label commonXYZ Node Fluc AveTemp CASE SNAP Para xyzD xyzU DataD DataU PredList PredCase PRED Method W
    % close all;
    clc;
    CASE = 10; % total cases
    Method = 0;  % using fluctuation (Method =1); otherwise, Method ~=1
    PredCase = PredList(PRED);
    AveTemp = cell2mat(arrayfun(@(k) mean(Common{k},2),1:CASE,'UniformOutput',false)); % %  Temporal Average % %
    [n,m]= size(Common{1});
    X = ones(n,m);
    % tic
    COVA = ones(m,m);
    for Lp=1:CASE
        if Lp ==PredCase
            continue
            tic
        else
            if Method == 1; %using fluctuation
                X = bsxfun(@minus,Common{Lp},AveTemp(:,Lp));
            else %using original data
                X = Common{Lp};
            end
            COVA = COVA.*((X'*X)./(n-1));
            tt = toc;
            clear X;
            fprintf('\n *********** Covariance for case %d is done with %5.2f sec *********** \n',Lp,tt);
        end
    end
    % toc
    
    % %====== Eigenvalues and eigenvectors ======
    [U,S] = eig(COVA);
    S = diag(S);
    [~, idx] = sort(S,'descend'); %sort by eigenvalue
    S = S(idx);
    U = U(:,idx);
    Energy = cumsum(S)./sum(S)*100;
    CutNUM = SNAP;
    %     CutNUM = length(Energy(Energy<=99.99));
    
    % % ====== k-Step Arnoldi Method for Eigenvalues and eigenvectors ======
    %     CutNUM = SNAP -10;
    fprintf('\n ***************************************************************** \n');
    fprintf('\n *********************** Using %d modes!! *********************** \n',CutNUM);
    fprintf('\n ***************************************************************** \n');
    
    %     [U,S] = eigs(COVA,CutNUM);
    %     S = diag(S);
    
    PODModes = cell(CASE,1);
    PODCoeff = cell(CASE,1);
    P = 0;
    
    for Lp = 1:CASE
        if Lp ==PredCase
            continue
            tic
        else
            if Method == 1; %using fluctuation
                if CutNUM >=600
                    for G=1:CutNUM
                        X(:,G) = Common{Lp}(:,G) - AveTemp(:,Lp);
                    end
                else
                    X = bsxfun(@minus,Common{Lp},AveTemp(:,Lp));
                end
            else %using original data
                X = Common{Lp};
            end
            PODModes{Lp} = X*U(:,1:CutNUM);
            PODCoeff{Lp} = U(:,1:CutNUM);
            
            L = zeros(1,CutNUM);
            for I = 1:CutNUM
                L(:,I) = norm(PODModes{Lp}(:,I));
            end
            if CutNUM>= 600
                for G=1:CutNUM
                    PODModes{Lp}(:,G) = PODModes{Lp}(:,G)./L(G);
                    PODCoeff{Lp}(:,G) = PODCoeff{Lp}(:,G).*L(G);
                end
            else
                PODModes{Lp} = bsxfun(@rdivide,PODModes{Lp},L);
                PODCoeff{Lp} = (bsxfun(@times,PODCoeff{Lp},L));
            end
            %     A = PODModes{Lp}*PODCoeff{Lp}';
            tt(Lp,1)= toc;
            fprintf('\n ****************Case %d is done with %5.2f sec************************** \n',Lp,tt(Lp,1));
        end
    end
    
    PODModes = PODModes(~cellfun('isempty',PODModes));
    PODCoeff = PODCoeff(~cellfun('isempty',PODCoeff));
    
    clc
    fprintf('\n ***************************************************************** \n');
    fprintf('\n ******************** Common-POD is complete! ******************** \n');
    fprintf('\n ***************************************************************** \n');
    
    
    %% ****** Emulation ******
    clearvars -except Common Label commonXYZ Node Fluc AveTemp CASE SNAP Para xyzD xyzU DataD DataU PredList PredCase PRED Method PODModes PODCoeff W
    %% =================== GP Process ================================== %%
    % % ********* Loading design matrix *********
    DES = ['..\Design Matrix_RD170.xlsx']; % design matrix file name
    filename  = fullfile(DES);  % combination of file name & folder dir
    DesignM = readtable(filename,'Range','A2:C8','ReadVariableNames',false);
    DesignM.Properties.VariableNames = {'Case','FuelRecess','Class'};
    DesignMatrix = table2array(DesignM);
    DesignMatrix(:,2) = DesignMatrix(:,2)/max(DesignMatrix(:,2));
    
    CutNUM = size(PODCoeff{1},2);
    
    Design = DesignMatrix(:,1:2);
    Xtest = Design(PredCase,:);
    Design(PredCase,:) = [];
    Xtrain = Design;
    
    
    ValiCase = Common{PredCase};
    nfile = SNAP;
    CASE0 = CASE-1;
    YtrainNEW = double(permute(reshape(cell2mat(PODCoeff),[SNAP,CASE0,CutNUM]),[2,3,1])); % training cases POD Coefficient
    
    
    YtestCoeff = zeros(nfile,CutNUM);
    for H = 1:nfile
        tic
        parfor (LP = 1:CutNUM,10)
            if (max(abs(YtrainNEW(:,LP,H)))>= 1e+3)
                YtrainNEW(:,LP,H) = YtrainNEW(:,LP,H).*(1e-3);
                gprMdl = fitrgp(Xtrain,YtrainNEW(:,LP,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern32',...
                    'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
                YtestCoeff(H,LP) = predict(gprMdl,Xtest);
                YtestCoeff(H,LP) = YtestCoeff(H,LP).*(1e+3);
            elseif (max(abs(YtrainNEW(:,LP,H)))>= 1e-3 && max(abs(YtrainNEW(:,LP,H)))<=1)
                YtrainNEW(:,LP,H) = YtrainNEW(:,LP,H).*(1e3);
                gprMdl = fitrgp(Xtrain,YtrainNEW(:,LP,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern32',...
                    'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
                YtestCoeff(H,LP) = predict(gprMdl,Xtest);
                YtestCoeff(H,LP) = YtestCoeff(H,LP).*(1e-3);
            elseif (max(abs(YtrainNEW(:,LP,H)))<= 1e-3)
                YtrainNEW(:,LP,H) = YtrainNEW(:,LP,H).*(1e6);
                gprMdl = fitrgp(Xtrain,YtrainNEW(:,LP,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern32',...
                    'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
                YtestCoeff(H,LP) = predict(gprMdl,Xtest);
                YtestCoeff(H,LP) = YtestCoeff(H,LP).*(1e-6);
            else
                gprMdl = fitrgp(Xtrain,YtrainNEW(:,LP,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern32',...
                    'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
                YtestCoeff(H,LP) = predict(gprMdl,Xtest);
            end
        end
        clear gprMdl
        fprintf('\n --------- GP for new %d snapshot coefficients is done--------------------- \n',H);
        toc
    end
    
    % % ****** POD Mode Weighting parameters ******
    
    W0 = zeros(CASE0,1);
    YtrainW = eye(CASE0);
    W1 = zeros(CASE0,1);
    %% ****** GPfit for new modes (KSPOD) ******
    for J=1:2
        tic
        if J==1
            parfor (H = 1:CASE0)
                gprMd2 =  fitrgp(Xtrain,YtrainW(:,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern52',...
                    'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
                W0(H,1) = predict(gprMd2,Xtest);
            end
        else
            tic
            parfor (H = 1:CASE0)
                gprMd2 =  fitrgp(Xtrain,W1(:,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern52',...
                    'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
                W0(H,1) = predict(gprMd2,Xtest);
            end
        end
        t = toc;
        fprintf('\n ****** Kernel Round %d: %4.6f sec ******\n\n',J,t);
        W1 = (W0)/sum((W0));
        W1 = diag(W1);
    end
    clc;
    W = diag(W1);
    
    
    %% ****** Reconstruction ******
    tic;
    for H=1:CASE0
        if H==1
            Prediction = bsxfun(@times,W(H), PODModes{H})*YtestCoeff';
        else
            Prediction = Prediction+bsxfun(@times,W(H), PODModes{H})*YtestCoeff';
        end
        fprintf('\n ****** Mode Prediciton for Case %d is done ******\n',H);
    end
    toc
    
    if Method == 1
        NewAveTemp = mean(ValiCase,2);
        Prediction = bsxfun(@plus,Prediction, NewAveTemp);
    end
    
    clc
    fprintf('\n ***************************************************************** \n');
    fprintf('\n ******************** Prediction is complete! ******************** \n');
    fprintf('\n ***************************************************************** \n');
    
    %% ===================== Save Prediction =====================
    clearvars -except Common Label commonXYZ Node Fluc AveTemp CASE SNAP Para xyzD xyzU DataD DataU PredList PredCase PRED Method PODModes PODCoeff Prediction ValiCase nfile W
    
    clc
    if Para == 7
        DirOut='../Output/Pressure/20180418';
        Tempfile  = ['D:\NewEmulationData\RD-170Group\PLTtemplate_P.plt'];
    elseif Para == 8
        DirOut='../Output/Temperature/20180418';
        Tempfile  = ['D:\NewEmulationData\RD-170Group\PLTtemplate_T.plt'];
    elseif Para == 9
        DirOut='../Output/Density/20180418';
        Tempfile  = ['D:\NewEmulationData\RD-170Group\PLTtemplate_D.plt'];
    else
        DirOut='../Output/Others';
    end
    %     fileFolder  = fullfile(DirOut,OutFile);
    eval(['OutFile=[''' 'Pred_Case=' num2str(PredCase) '' '''];']);
    fileFolder  = fullfile(DirOut,OutFile);
    [s,mess,messid] = mkdir(fileFolder);
    
    if max(max(commonXYZ))<= 10
        commonXYZ = commonXYZ.*1000; %unit mm
    else
        commonXYZ = commonXYZ;
    end
    
        
    %% ========= Making Template for PLT =========
    eval(['OutFile=[''' 'Pred_Case=' num2str(PredCase) '' '''];']);
    fileDEST  = fullfile(DirOut,OutFile);
    [s,mess,messid] = mkdir(fileDEST);
    
    fprintf('****** Save Emulation ****** \n\n');
    DATA3 = Node{3};
    if max(max(commonXYZ))<= 1
        commonXYZ = commonXYZ.*1000; %unit mm
    else
        commonXYZ = commonXYZ;
    end
    parfor (Snapshot = 1:nfile)
        tic
        filename = ['Pred_Case=' num2str(PredCase) '_' num2str(Snapshot) '.dat'];
        fileDEST  = fullfile(DirOut,OutFile,filename);
        fileID = fopen(fileDEST,'wt');
        
        fprintf(fileID,'Title = "Reconstruction" \n');
        fprintf(fileID,'VARIABLES = "X" \n"Y" \n');
        fprintf(fileID,'"Simulation" \n');
        fprintf(fileID,'"Emulation" \n');
        
        
        fprintf(fileID,'ZONE T="Triangulation"\n');
        fprintf(fileID,'STRANDID=0, SOLUTIONTIME=0 \n');
        fprintf(fileID,'Nodes=75327, Elements=149083, ZONETYPE=FETriangle\n');
        fprintf(fileID,'DATAPACKING=POINT\n');
        fprintf(fileID,' DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n');
        
        for ijk = 1:size(commonXYZ,1)
            fprintf(fileID,'%e\t',commonXYZ(ijk,:),ValiCase(ijk,Snapshot),Prediction(ijk,Snapshot));
            fprintf(fileID,'\n');
        end
        
        
        
        for ijk  =1:size(DATA3,1)
            fprintf(fileID,'%d\t',DATA3(ijk,:));
            fprintf(fileID,'\n');
        end
        
        fclose(fileID);
        toc
    end
    
    clearvars -except Common Label commonXYZ Node CASE SNAP Para xyzD xyzU DataD DataU PredList PredCase PRED Method PODModes PODCoeff Prediction ValiCase nfile W
    
    
end

