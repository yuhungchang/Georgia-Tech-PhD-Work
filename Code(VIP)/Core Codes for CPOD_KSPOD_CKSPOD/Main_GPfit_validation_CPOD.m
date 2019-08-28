%% =================== GP Process (CPOD) ================================== %%
clc;
SNAP = Info.SNAP;
CASE = Info.SNAP;
% % ********* Loading design matrix *********
filename = ['Design Matrix_new 30 cases.xlsx']; % design matrix file name 
DesignM = readtable(filename,'Range','B2:I70','ReadVariableNames',false);
DesignM.Properties.VariableNames = {'h','theta','dL','K','Cluster','h_real','theta_real','dL_real'};

% PredList =[67 66 64];  % prediction case (for validation)
PredList =[61 62 64 65 66 67 68 69 63];  % prediction case (for validation)

% nfile = SNAP; % snapshot number for prediction

clear DATA;
clc;
for LL= 6:8 %7:length(PredList)
    close all; clc
    clearvars -except Method Para List DATA XYZ PODCoeff PODModes CutNUM PredList LL Info W DesignM Core nfile
        
    %% ====== Loading Validation Case ======
    PredCase = PredList(LL);
    fprintf('\n *********************************************** \n');
    fprintf('\n ********* Loading Validation Case %d ********** \n',PredCase);
    fprintf('\n *********************************************** \n');
    
    
    clear ValiCase
    List = Info.List; Para = Info.Para; Var= Info.Var; Err = Info.Err; Kt =Info.Kt; I=Info.I; J=Info.J;
    SKIP = Info.SKIP; SNAP = Info.SNAP; CASE = Info.CASE; d1 = Info.d1; d2=Info.d2;
    
    k1 = int2str(PredList(LL));
%     ValiCase = zeros(size(DATA{1}));
    tic
    Folder =['../case',k1,'_plt'];% file name combination
    [ValiCase,~] = LoadingPLT(Kt, I,J , Var,Err,SKIP,Para,d1,d2,Folder);
    ValiCase = double(ValiCase);
    clear XYZorg
    toc
    
    fprintf('\n --------- **** Finish loading test case %s for validation **** --------------------- \n',int2str(PredList(LL)));
    
    %     clc;
    clearvars -except Method Para List DATA XYZ PODCoeff PODModes CutNUM PredList LL Info W DesignM Core ValiCase nfile
    %% ====== Prepare for GPfit ======
   
    PredCase = PredList(LL);
    % % ********* Choosing training cases and testing cases*********
    SNAP = Info.SNAP; CASE = Info.CASE; SKIP = Info.SKIP;
    Xtrain = table2array(DesignM(List,[1:3,5])); % training cases
    Xtest = table2array(DesignM(PredCase,[1:3,5])); % testing (or prediction) case
    YtrainNEW = double(permute(reshape(cell2mat(PODCoeff),[SNAP,CASE,CutNUM]),[2,3,1])); % training cases POD Coefficient
    %% =================== GP with POD ================================== %%
    Snapname = (1:(SKIP+1):1000)';
    ST = 1;
    TimeSTEP = 20; %SNAP;
    nfile = TimeSTEP+ST-1;     
    if nfile>length(Snapname)
        nfile = length(Snapname);
        ST= nfile - TimeSTEP+1;
    end
    
    YtestCoeff = zeros(TimeSTEP,CutNUM);
    % ******************** Coefficient GP process ********************
    for H = ST:nfile
        tic
        parfor (LP = 1:CutNUM,Core)
            gprMdl = fitrgp(Xtrain,YtrainNEW(:,LP,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern32',...
                'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
            YtestCoeff(H,LP) = predict(gprMdl,Xtest);
            
            % ========== Use the part below if Max and Min of YtrainNEW is huge ==========
            %             if (max(abs(YtrainNEW(:,LP,H)))>= 1e+3)
            %                 YtrainNEW(:,LP,H) = YtrainNEW(:,LP,H).*(1e-4);
            %                 gprMdl = fitrgp(Xtrain,YtrainNEW(:,LP,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern32',...
            %                     'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
            %                 YtestCoeff(H,LP) = predict(gprMdl,Xtest);
            %                 YtestCoeff(H,LP) = YtestCoeff(H,LP).*(1e+4);
            %             elseif(max(abs(YtrainNEW(:,LP,H)))<= 1e-3)
            %                 YtrainNEW(:,LP,H) = YtrainNEW(:,LP,H).*(1e+4);
            %                 gprMdl = fitrgp(Xtrain,YtrainNEW(:,LP,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern32',...
            %                     'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
            %                 YtestCoeff(H,LP) = predict(gprMdl,Xtest);
            %                 YtestCoeff(H,LP) = YtestCoeff(H,LP).*(1e-4);
            %             else
            %                 gprMdl = fitrgp(Xtrain,YtrainNEW(:,LP,H),'BasisFunction','pureQuadratic','FitMethod','fic','KernelFunction','ardmatern32',...
            %                     'PredictMethod','fic','Sigma',1e-15,'Standardize',1,'DistanceMethod','accurate','Optimizer','fmincon','SigmaLowerBound',1e-9);
            %                 YtestCoeff(H,LP) = predict(gprMdl,Xtest);
            %             end
            % ============================================================================
            
        end
        clear gprMdl
        fprintf('\n --------- GP for new %d snapshot coefficients is done--------------------- \n',Snapname(H));
        toc
    end
    
    
    %% ****** Reconstruction ****** 
    tic;
    Prediction = PODModes*YtestCoeff';
 
    toc
    if Method == 1 % using fluctuation  (if Method ==1, it means the prediction only carries fluctuation information)
        NewAveTemp = mean(ValiCase,2);
        if CutNUM <= 600
            Prediction = bsxfun(@plus,Prediction, NewAveTemp); % adding time average information back
        else
            for ss = 1:SNAP
                Prediction(:,ss) = Prediction(:,ss) + NewAveTemp; % adding time average information back
            end
        end
    end
    clc
    fprintf('\n ***************************************************************** \n');
    fprintf('\n ******************** Prediction is complete! ******************** \n');
    fprintf('\n ***************************************************************** \n');
    
    
    %% ===================== Save Prediction ===================== %%
    clc;
    if Para ==7
        DirOut=['../OutPut/Pressure/CPOD_20180413_skip=' num2str(SKIP) '']; 
        ValiCase = ValiCase.*(1e-4);
        Prediction = Prediction.*(1e-4);
        Tempfile = fopen('PLTTemplate_Pressure.plt');
    elseif Para == 8
        DirOut=['../OutPut/Temperature/CPOD_20180413_skip=' num2str(SKIP) ''];         
        Tempfile = fopen('PLTTemplate_Temperature.plt');
    elseif Para == 9
        DirOut=['../OutPut/Density/CPOD_20180413_skip=' num2str(SKIP) ''];
        Tempfile = fopen('PLTTemplate_Density.plt');
    else         
        DirOut=['../OutPut/Others/CPOD_20180413_skip=' num2str(SKIP) ''];
        Tempfile = fopen('PLTTemplate.plt');
    end
    eval(['OutFile=[''' 'CPOD_case' num2str(PredCase) '_mode=' num2str(CutNUM) '' '''];']); 
    
    fileFolder  = fullfile(DirOut,OutFile);
    [~,mess,messid] = mkdir(fileFolder); clc;
    
    %% **** Loading Header and other information for binary file header ****    
    
    A = fread(Tempfile,'*single');
    B = flipud(A);
    idx = find(B==299); % '299' = zone
    Var0 = 5; I=33; J=33; IJ = I*J;
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
    Snapname = (1:(SKIP+1):1000)';
    tic
    for H = ST:nfile
        RawData = [XYZ ValiCase(:,H) Prediction(:,H)];
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
        Snapshot = Snapname(H);
        pltDEST = ['Case_' num2str(PredCase) '_pred_' num2str(Snapshot) '.plt'];
        if exist(pltDEST, 'file')==2
            delete(pltDEST);
        end
        eval(['!rename' 32 'raw.bin' 32 pltDEST]);
        movefile(pltDEST,fileFolder);
        clear DData3 DData2 DData1 binDEST pltDEST
        
    end
    toc
    fprintf('\n ****************************************************************\n');
    fprintf('\n ********** Prediction data writing for case %d is done ********* \n',PredCase);
    fprintf('\n ****************************************************************\n');
    % delete(gcp); 
    clearvars -except Method Para List DATA XYZ PODCoeff PODModes CutNUM PredList LL Info W DesignM Core nfile
        
end


