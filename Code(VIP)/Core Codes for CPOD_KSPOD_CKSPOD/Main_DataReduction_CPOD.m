%% ===================== Data Reduction by Proper Orthogonal Decomposition (POD) =====================
clearvars -except DATA XYZ List Para SKIP CASE SNAP CutNUM Start Pt Core Info TestCase XYZorg PredList W Core
clc;
CASE = Info.CASE;  SNAP = Info.SNAP;
PODCoeff0 = cell(CASE,1); PODModes = cell(CASE,1); Percent = cell(CASE,1);


fprintf('\n ************************************************************ \n');
fprintf('\n ********* Creating Common Covariance over %d cases ********* \n', CASE);
fprintf('\n ************************************************************ \n');
Method = 0;  % using fluctuation (Method =1); otherwise, Method ~=1

AveTemp = cell2mat(arrayfun(@(k) mean(DATA{k},2),1:CASE,'UniformOutput',false)); % %  Temporal Average % %
[n,m]= size(DATA{1});
X = []; %ones(n,m);
% tic
% COVA = ones(SNAP,SNAP);
for Lp=1:CASE
    tic
    if Method == 1; %using fluctuation
        X0 = bsxfun(@minus,DATA{Lp},AveTemp(:,Lp));
    else %using original data
        X0 = DATA{Lp};
    end
    X = cat(2,X,X0);
    tt = toc;
    clear X0;
    fprintf('\n *********** Covariance for case %d is done with %5.2f sec *********** \n',Lp,tt);
end

tic
COVA = ((X'*X)./(n-1));
t1=toc;
 
tic
[~,S] = eig(COVA);
S = diag(S);
[~, idx] = sort(-1*S); %sort by eigenvalue
S = S(idx);
% U = U(:,idx);
Energy = cumsum(S)/sum(S)*100;
CutNUM = length(Energy(Energy<=99.5));
t2=toc;
% CutNUM = 800;
fprintf('\n ************************************************************ \n');
fprintf('\n ********* SVD and cut-off at 99.5%% with %5.2f sec ********* \n',(t1+t2));
fprintf('\n ************************************************************ \n');

clear U S
% CutNUM = SNAP-1;
[U,S] = eigs(COVA,CutNUM);
S = diag(S);

% PODenery(:,Lp)=Energy;

% % *********** calculating modes ******** % %
fprintf('Calculate POD modes & Coefficients ...\n');
clear PODCoeff0 PODModes
tic
PODModes= X*U(:,1:CutNUM);

PODCoeff0 = PODModes'*X;
if size(X,2)>= 500
    for J=1:size(PODCoeff0,2)
        PODCoeff0(:,J) = PODCoeff0(:,J)./(S(1:CutNUM)*(n-1));
    end
    PODCoeff0 = PODCoeff0';
else
    PODCoeff0 = (bsxfun(@rdivide,PODCoeff0,(S(1:CutNUM)*(n-1))))';
end
tt= toc;
fprintf('\n ****************Case %d is done with %5.2f sec************************** \n',Lp,tt);
clear X S U;

CASE = Info.CASE;
PODCoeff = mat2cell(PODCoeff0,SNAP*ones(1,CASE),CutNUM);
clear PODCoeff0
% PODModes = PODModes(~cellfun('isempty',PODModes));
% PODCoeff0 = PODCoeff0(~cellfun('isempty',PODCoeff0));
 
% A = PODModes*PODCoeff0';
clc;
fprintf('\n ***************************************************************** \n');
fprintf('\n ******************** Common-KSPOD is complete! ****************** \n');
fprintf('\n ***************************************************************** \n');

