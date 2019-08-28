function [XYZDown,XYZUp,XYZR,XYZUpB,Downstream,Upstream,Recess,UpstreamB] = PLTloaderRD170(I,J,Var,Para,Range,fileDEST)
% %  % ************************************  %  %
% %  % *********  By Yu-Hung Chang  *******  %  %
% %  % ***********  2016-Aug-29  **********  %  %
% %  % ************************************  %  %

IJ = I*J;
d1 = Range(1);
d2 = Range(2);
d3 = Range(3);
r1 = Range(4); %Length of
r2 = Range(5);
r3 = Range(6);
err = 0.00005;

% Ln = 0.1035;
% %  % **** Err: lines between zones ******  %  %
% %  % **** Zone: Zones in one CFD case ***  %  %
% %  % **** Var: number of variables  *****  %  %
% %  % **** countNUM: number of snapshots *  %  %
% %  % **** CASE: number of Cases *********  %  %


fileID = fopen(fileDEST);
A = fread(fileID,'*single')';
fclose(fileID);
A = A';
[Ind,~] = find(A==299);
% L0 = diff(Ind);
L1 = max(diff(Ind));
[Zone,~] = find(diff(Ind)==L1);
Zone = Zone(1)-1;
B = flipud(A);

% B = flipud(A'); clear A;


num = L1*Zone;
B(num+1:end,:) =[];

RawData = reshape(flipud(B),[],Zone); 
L2 = L1-IJ*Var;

RawData(1:L2,:) =[];
if RawData(2,1)==0
    OutPut = (reshape(permute(reshape(RawData',Zone,Var,IJ),[2 3 1]),Var,[]))';
else
    OutPut = reshape(permute(reshape(RawData',Zone,IJ,Var),[2 1 3]),[],Var);
end
% clear DATA A B
OutPut(:,3)=1;
OutPut((OutPut(:,1)>=(r1+err)),3) = NaN; %downstream = NaN

%% ========= Upstream =========
OutPut1 = permute(reshape(OutPut,[IJ,Zone,Var]),[1,3,2]);
idxToRemove = any(any(isnan(OutPut1),2),1);
OutPut1(:,:,idxToRemove) = []; %upstream
OutPut1 = reshape(permute(OutPut1,[2 1 3]),Var,[])';  %upstream

OutPut1(:,3)=1;
OutPut1((OutPut1(:,1)< d2-err*1.2),3) = NaN;
% OutPut1((OutPut1(:,2)> r2+0.0005),3) = NaN;
NewZone = size(OutPut1,1)/IJ;
OutPut1 = permute(reshape(OutPut1,[IJ,NewZone,Var]),[1,3,2]);
idxToRemove = any(any(isnan(OutPut1),2),1);
OutPut1(:,:,idxToRemove) = []; %upstream
OutPut1 = reshape(permute(OutPut1,[2 1 3]),Var,[])';  %upstream overall

OutPut2 = OutPut1; 
OutPut2(:,3)=1;
OutPut2((OutPut2(:,2)< r2-err*1.2),3) = NaN; % below recess area
NewZone = size(OutPut2,1)/IJ; 
OutPut2 = permute(reshape(OutPut2,[IJ,NewZone,Var]),[1,3,2]);
idxToRemove = any(any(isnan(OutPut2),2),1);
OutPut2(:,:,idxToRemove) = [];  
OutPut2 = reshape(permute(OutPut2,[2 1 3]),Var,[])';  %upstream below recess area
Upstream = OutPut2(:,Para);
XYZUp = OutPut2(:,1:2); clear  OutPut2 


OutPut2 = OutPut1; 
OutPut2(:,3)=1;
OutPut2((OutPut2(:,2)< r2+err*1.2),3) = NaN; % below recess area
OutPut2((OutPut2(:,2)> r3+err*1.2),3) = NaN; % below recess area
NewZone = size(OutPut2,1)/IJ; 
OutPut2 = permute(reshape(OutPut2,[IJ,NewZone,Var]),[1,3,2]);
idxToRemove = any(any(isnan(OutPut2),2),1);
OutPut2(:,:,idxToRemove) = [];  
OutPut2 = reshape(permute(OutPut2,[2 1 3]),Var,[])';  %upstream below recess area
Recess = OutPut2(:,Para);
XYZR = OutPut2(:,1:2); clear  OutPut2 

OutPut2 = OutPut1; 
OutPut2(:,3)=1; 
OutPut2((OutPut2(:,2)< r3+err*1.2),3) = NaN; % below recess area
NewZone = size(OutPut2,1)/IJ; 
OutPut2 = permute(reshape(OutPut2,[IJ,NewZone,Var]),[1,3,2]);
idxToRemove = any(any(isnan(OutPut2),2),1);
OutPut2(:,:,idxToRemove) = [];  
OutPut2 = reshape(permute(OutPut2,[2 1 3]),Var,[])';  %upstream below recess area
UpstreamB = OutPut2(:,Para);
XYZUpB = OutPut2(:,1:2); clear  OutPut2 


% ======== separate the recess part =======
% Upstream0 = OutPut1(:,Para);
% XYZUp0 = OutPut1(:,1:2);
% 
% XYZUp0(:,3) = ones(length(XYZUp0),1);
% XYZUp0((XYZUp0(:,2)< r2),3) = NaN;
% 
% Upstream = Upstream0(any(isnan(XYZUp0),2),:); 
% XYZUp = XYZUp0(any(isnan(XYZUp0),2),:); 
% XYZUp(:,3)=[];
% 
% Recess0 = Upstream0(~any(isnan(XYZUp0),2),:); 
% XYZR0 = XYZUp0(~any(isnan(XYZUp0),2),:); 
% 
% XYZR0(:,3)= ones(length(XYZR0),1);
% XYZR0((XYZR0(:,2)< r3),3) = NaN;
% 
% Recess = Recess0(any(isnan(XYZR0),2),:); 
% XYZR = XYZR0(any(isnan(XYZR0),2),:); 
% XYZR(:,3)=[];
% 
% UpstreamB = Recess0(~any(isnan(XYZR0),2),:); 
% XYZUpB = XYZR0(~any(isnan(XYZR0),2),:); 


%% ========= Downstream =========
OutPut2 = permute(reshape(OutPut,[IJ,Zone,Var]),[1,3,2]);
idxToRemove = any(any(isnan(OutPut2),2),1);
OutPut2(:,:,~idxToRemove) = []; %downstream
OutPut2 = reshape(permute(OutPut2,[2 1 3]),Var,[])';  %downstream

OutPut2(:,3)=1;
OutPut2((OutPut2(:,1)> d1+err)|(OutPut2(:,2)> d3+err),3) = NaN;
NewZone = size(OutPut2,1)/IJ;
OutPut2 = permute(reshape(OutPut2,[IJ,NewZone,Var]),[1,3,2]);
idxToRemove = any(any(isnan(OutPut2),2),1);
OutPut2(:,:,idxToRemove) = []; %downstream
OutPut2 = reshape(permute(OutPut2,[2 1 3]),Var,[])';  %downstream

Downstream = OutPut2(:,Para);
XYZDown = OutPut2(:,1:2);


% M1 = [XYZDown Downstream]; 
% M2 = [XYZUp Upstream];

end
