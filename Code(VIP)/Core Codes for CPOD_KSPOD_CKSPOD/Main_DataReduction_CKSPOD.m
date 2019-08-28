%% ===================== Data Reduction by Proper Orthogonal Decomposition (POD) =====================
% clearvars -except DATA XYZ List Para SKIP CASE SNAP CutNUM Start Pt Core Info TestCase XYZorg PredList W Core
clc;
CASE = Info.CASE;  SNAP = Info.SNAP;
Percent = cell(CASE,1);

fprintf('\n ************************************************************ \n');
fprintf('\n ********* Creating Common Covariance over %d cases ********* \n', CASE);
fprintf('\n ************************************************************ \n');

Method = 0;  % using fluctuation (Method =1); otherwise, Method ~=1
SafeFile = 0; % choose "SafeFile = 1" to save CKSPOD information (the first 10 POD modes only)

AveTemp = cell2mat(arrayfun(@(k) mean(DATA{k},2),1:CASE,'UniformOutput',false)); % %  Temporal Average % %
[n,m]= size(DATA{1});
X = ones(n,m);
% tic
COVA = ones(SNAP,SNAP);
for G=1:CASE
    tic
    if Method == 1; %using fluctuation
        X = bsxfun(@minus,DATA{G},AveTemp(:,G));
    else %using original data
        X = DATA{G};
    end
    COVA = COVA.*((X'*X)./(n-1));
    tt = toc;
    clear X;
    fprintf('\n *********** Covariance for case %d is done with %5.2f sec *********** \n',G,tt);
end
% toc

% CutNUM = 800;
fprintf('\n ********************************************* \n');
fprintf('\n ********* SVD and cut-off at 99.9%% ********* \n');
fprintf('\n ********************************************* \n');

[~,S] = eig(COVA);
S = diag(S);
[~, idx] = sort(-1*S); %sort by eigenvalue
S = S(idx); % U = U(:,idx);
EnergyCKS = cumsum(S)/sum(S)*100;
% CutNUM = length(EnergyCKS(EnergyCKS<99.5));
clear U S;
% CutNUM = 150;
CutNUM = SNAP-1;
[U,S] = eigs(COVA,CutNUM);
S = diag(S);

PODModes = cell(CASE,1);
PODCoeff = cell(CASE,1);
fprintf('\n ***************************************************** \n');
fprintf('\n ********* Creating POD modes & coefficients ********* \n');
fprintf('\n ***************************************************** \n');
tt = zeros(CASE,1);
% parfor (Lp =1:CASE,Core)
for G =1:CASE
    tic
    if Method == 1; %using fluctuation
        if CutNUM >=600
            for G=1:SNAP
                X(:,G) = DATA{G}(:,G) - AveTemp(:,G);
            end
        else
            X = bsxfun(@minus,DATA{G},AveTemp(:,G));
        end
    else %using original data
        X = DATA{G};
    end
    PODModes{G} = X*U(:,1:CutNUM);
    PODCoeff{G} = U(:,1:CutNUM);
    
    L = zeros(1,CutNUM);
    for I = 1:CutNUM
        L(:,I) = norm(PODModes{G}(:,I));
    end
    if CutNUM>= 600
        for G=1:CutNUM
            PODModes{G}(:,G) = PODModes{G}(:,G)./L(G);
            PODCoeff{G}(:,G) = PODCoeff{G}(:,G).*L(G);
        end
    else
        PODModes{G} = bsxfun(@rdivide,PODModes{G},L);
        PODCoeff{G} = bsxfun(@times,PODCoeff{G},L);
    end
    A = PODModes{G}*PODCoeff{G}';
    tt(G,1)= toc;
    fprintf('\n ****************Case %d is done with %5.2f sec************************** \n',G,tt(G,1));
    
end

clear X L;
PODModes = PODModes(~cellfun('isempty',PODModes));
PODCoeff = PODCoeff(~cellfun('isempty',PODCoeff));

% Lp = 10;
% A = PODModes{Lp}*PODCoeff{Lp}';
clc;
fprintf('\n ***************************************************************** \n');
fprintf('\n ******************** Common-KSPOD is complete! ****************** \n');
fprintf('\n ***************************************************************** \n');

%% ===================== Save CKSPOD POD modes ===================== %%
% clearvars -except Method Normal CutOff Para List DATA XYZ PODCoeff PODModes CutNUM Info Core EnergyCKS SafeFile 
SKIP= Info.SKIP; SNAP = Info.SNAP; CASE = Info.CASE;
if SafeFile == 1 % choose to save CKSPOD information (the first 10 POD modes only)
    clc;
    if Para ==7
        DirOut=['../OutPut/POD/Pressure/CKSPOD_skip=' num2str(SKIP) ''];
    elseif Para == 8
        DirOut=['../OutPut/POD/Temperature/CKSPOD_skip=' num2str(SKIP) ''];
    elseif Para == 9
        DirOut=['../OutPut/POD/Density/CKSPOD_skip=' num2str(SKIP) ''];
    else
        DirOut=['../OutPut/POD/Others/CKSPOD_skip=' num2str(SKIP) ''];
    end
    [~,mess,messid] = mkdir(DirOut); clc;
    CutOff = 10; Normal =1;
    for G = 1:CASE
        if Method ==1
            filename1 = ['Fluc_MODE_Case=' num2str(G) '.plt'];
            filename2 = ['Fluc_Coeff_Case=' num2str(G) '.dat'];
        else
            filename1 = ['All_MODE_Case=' num2str(G) '.plt'];
            filename2 = ['All_Coeff_Case=' num2str(G) '.dat'];
        end
        fileDEST1  = fullfile(DirOut,filename1);
        fileDEST2  = fullfile(DirOut,filename2);
        
        WritingPODplt(CutOff,G,SNAP,XYZ,PODModes{G},PODCoeff{G},fileDEST1,fileDEST2,Normal);
        %         WritingPOD(CutOff,num,IJK,XYZ,MODES,Coeff,fileDEST1,fileDEST2,CaseNUM);
    end
    %     clear MODES Coeff
end

%%  ================ Save POD information for raw data ========================
% clearvars -except Method Normal CutOff Para List DATA XYZ PODCoeff PODModes CutNUM Info Core EnergyCKS SafeFile
SKIP= Info.SKIP; SNAP = Info.SNAP; CASE = Info.CASE;
EnergyRaw = zeros(SNAP,CASE);
T1 = tic;
if SafeFile == 1 % choose to save CKSPOD information (the first 10 POD modes only)
    clc;
    if Para ==7
        DirOut=['../OutPut/POD/Pressure/POD_skip=' num2str(SKIP) ''];
    elseif Para == 8
        DirOut=['../OutPut/POD/Temperature/POD_skip=' num2str(SKIP) ''];
    elseif Para == 9
        DirOut=['../OutPut/POD/Density/POD_skip=' num2str(SKIP) ''];
    else
        DirOut=['../OutPut/POD/Others/POD_skip=' num2str(SKIP) ''];
    end
    [~,mess,messid] = mkdir(DirOut); clc;
    CutOff = 10; Normal =1;     
    
    for G = 1:CASE
        RawData = DATA{G};
        tic
        if Method == 1; %using fluctuation
            AveTemp = mean(RawData,2);
            X = bsxfun(@minus,RawData,AveTemp);
        else %using original data
            X = RawData;
        end
        clear RawData
        X = double(X);
        %     X = bsxfun(@minus,X,mean(X,2));
        [n,~]= size(X);
        COVA = ((X'*X)./(n-1));
        
        [~,S] = eig(COVA);
        S = diag(S);
        [~, idx] = sort(-1*S); %sort by eigenvalue
        S = S(idx); % U = U(:,idx);
        EnergyRaw(:,G) = cumsum(S)/sum(S)*100;
        EE = EnergyRaw(:,G) ;
        CutNUM = length(EE(EE<99)); clear EE
        clear U S;
        % CutNUM = 150;     % CutNUM = SNAP-100;
        [U,S] = eigs(COVA,CutNUM);
        S = diag(S);
        
        % % *********** calculating modes ******** % %
        fprintf('Calculate POD modes & Coefficients ...\n');
        tic
        MODES = X*U(:,1:CutNUM);
        
        COEFF = MODES'*X;
        if size(X,2)>= 500
            for J=1:size(COEFF,2)
                COEFF(:,J) = COEFF(:,J)./(S(1:CutNUM)*(n-1));
            end
            COEFF = COEFF';
        else
            COEFF = (bsxfun(@rdivide,COEFF,(S(1:CutNUM)*(n-1))))';
        end
        tt=  toc;
        fprintf('\n **************** POD for Case %d is done with %5.2f sec************************** \n',G,tt);
        Normal =1;
        
        if Method ==1
            filename1 = ['Fluc_MODE_Case=' num2str(G) '.plt'];
            filename2 = ['Fluc_Coeff_Case=' num2str(G) '.dat'];
        else
            filename1 = ['All_MODE_Case=' num2str(G) '.plt'];
            filename2 = ['All_Coeff_Case=' num2str(G) '.dat'];
        end
        fileDEST1  = fullfile(DirOut,filename1);
        fileDEST2  = fullfile(DirOut,filename2);
        
        WritingPODplt(CutOff,G,SNAP,XYZ,MODES,COEFF,fileDEST1,fileDEST2,Normal); 
        
    end
    
end
% T2 = tic; 
% Min = floor((T1-T2)/60); Sec = (T1-T2)-Min*60; 
clc
% fprintf('\n --------- Total time of saving POD files for %d cases: %d min %5.2f sec.  -------------------- \n',CASE,Min,Sec);
% fprintf('\n --------- Total time of data reduction process: %d min %5.2f sec.  -------------------- \n',Min,Sec);

% clearvars -except Method Normal CutOff Para List DATA XYZ PODCoeff PODModes CutNUM Info Core EnergyCKS EnergyRaw SafeFile



