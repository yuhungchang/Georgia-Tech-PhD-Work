function [PODcoef,PODmode,Energy,U,S] = POD(Data0,CutNUM,Volume,COND)
%  ************************************** %
%  ***********  June-03, 2017 ************ %
%  ********** By Yu-Hung Chang ********** %
%  ************************************** % 
 
close all; clc; 

[n,m]=size(Data0);


if nargin <4 
    COND = 1;  %default setting: Data are "fluctuations"
elseif nargin < 3     % If there's no volume size information
    Volume = ones(n,1);   
elseif nargin <2
    CutNUM = m;
end

% The cut-off number should be smaller or equal to the number of snapshots
if CutNUM > m
    CutNUM = m;
end

if COND==1 %True: Data0 (input) is fluctuations
    Data = Data0;
    clear Data0
else %False: Data0(input) is NOT fluctuations
    %% ===================== Spatial and Temporal Averages  ===================== %%
    clc;
    fprintf('\n --------- Calculating Spatial and Temporal Averages--------------------- \n');
    AveTemp = mean(Data0,2); % %  Temporal Average % %     
    Data = bsxfun(@minus,Data0,AveTemp); 
    fprintf('\n --------- Spatial and Temporal Average Calculation Is Done!!!--------------------- \n');    
    clear Data0
end 

X = bsxfun(@times,Data,Volume);
X = double(X);


fprintf('\n ********** Calculating Covariance or Auto-Correlation Matrix ********** \n');
tic;
COVA = X'*X./(n-1);
fprintf('\n ********** Covariance calculation is done ********** \n'); 

if CutNUM < m  % Arnoldi Approximation
    [U,S] = eigs(COVA,CutNUM);
    S = diag(S);

else   % Full Rank 
    [U,S] = eig(COVA);
    S = diag(S);
    [~, idx] = sort(-1*S); %sort by eigenvalue
    S = S(idx);
    U = U(:,idx); 
 
end
t0 = toc;
fprintf('\n****** SVD calculation processing time = %g sec ******\n',t0);
Energy = cumsum(S)/sum(S)*100;

% % *********** calculating modes ******** % %
fprintf('Calculate POD modes & Coefficients ...\n');
tic
PODmode = Data*U(:,1:CutNUM);
t1 = toc;
fprintf('\n****** POD mode calculation processing time = %g sec ******\n',t1);
 

tic;
if size(Data,2)>= 500
    
    for J=1:size(Data,2)
        X0(:,J) = Data(:,J).*Volume.^2;
    end
    
else
    X0 =bsxfun(@times,Data,Volume.^2);
end  

PODcoef = PODmode'*X0; 
if size(Data,2)>= 500
    for J=1:size(PODcoef,2)
        PODcoef (:,J) = PODcoef(:,J)./(S(1:CutNUM)*(n-1));
    end
    PODcoef = PODcoef';
else
    PODcoef = (bsxfun(@rdivide,PODcoef,(S(1:CutNUM)*(n-1))))';
end


t2 = toc;
fprintf('\n****** POD coefficient processing time = %g sec ******\n',t2);

% %     Recon = (PODcoef *PODmode')';
% %     Error = norm((Data -Recon),'fro')/norm(Data ,'fro');
 

end