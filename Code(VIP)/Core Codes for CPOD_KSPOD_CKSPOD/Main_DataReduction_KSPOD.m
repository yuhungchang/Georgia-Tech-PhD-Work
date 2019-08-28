%% ===================== Data Reduction by Proper Orthogonal Decomposition (POD) =====================
% clearvars -except DATA XYZ List Para SKIP CASE SNAP CutNUM Start Pt Core Info TestCase XYZorg PredList W Core
clc;
CASE = Info.CASE;  SNAP = Info.SNAP;
PODCoeff = cell(CASE,1); PODModes = cell(CASE,1); Percent = cell(CASE,1);


fprintf('\n ************************************************************ \n');
fprintf('\n ********* Creating Covariance over %d cases ********* \n', CASE);
fprintf('\n ************************************************************ \n');
Method = 0;  % using fluctuation (Method =1); otherwise, Method ~=1

AveTemp = cell2mat(arrayfun(@(k) mean(DATA{k},2),1:CASE,'UniformOutput',false)); % %  Temporal Average % %
[n,m]= size(DATA{1});
X = ones(n,m);
% tic

PODModes = cell(CASE,1);
PODCoeff = cell(CASE,1);
CutNUM  = SNAP-2;
PODenery = zeros(CutNUM,CASE);

for Lp=1:CASE
    tic
    if Method == 1; %using fluctuation
        X = bsxfun(@minus,DATA{Lp},AveTemp(:,Lp));
    else %using original data
        X = DATA{Lp};
    end
    X = double(X); 
%     X = bsxfun(@minus,X,mean(X,2));
    COVA = ((X'*X)./(n-1)); 
    
    
    [U,S] = eigs(COVA,CutNUM);
    S = diag(S);      
    Energy = cumsum(S)/sum(S)*100;
    PODenery(:,Lp)=Energy;
    
    % % *********** calculating modes ******** % %
    fprintf('Calculate POD modes & Coefficients ...\n');
    tic
    PODModes{Lp} = X*U(:,1:CutNUM); 
     
    PODCoeff{Lp} = PODModes{Lp}'*X;
    if size(X,2)>= 500
        for J=1:size(PODCoeff{Lp},2)
            PODCoeff{Lp} (:,J) = PODCoeff{Lp}(:,J)./(S(1:CutNUM)*(n-1));
        end
        PODCoeff{Lp} = PODCoeff{Lp}';
    else
        PODCoeff{Lp} = (bsxfun(@rdivide,PODCoeff{Lp},(S(1:CutNUM)*(n-1))))';
    end 
    tt(Lp,1)= toc;
    fprintf('\n ****************Case %d is done with %5.2f sec************************** \n',Lp,tt(Lp,1));
     
end
 
PODModes = PODModes(~cellfun('isempty',PODModes));
PODCoeff = PODCoeff(~cellfun('isempty',PODCoeff));
% clearvars -except PODenery DATA XYZ List Para SKIP CASE SNAP CutNUM Start Pt Core Info TestCase XYZorg PredList W Core PODModes PODCoeff Method

% Lp = 10;
% A = PODModes{Lp}*PODCoeff{Lp}';
clc;
fprintf('\n *********************************************************************** \n');
fprintf('\n ******************** Kernel-smoothedPOD is complete! ****************** \n');
fprintf('\n *********************************************************************** \n');

