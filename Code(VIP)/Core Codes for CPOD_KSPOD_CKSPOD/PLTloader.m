function [Results,XYZ] = PLTloader(Zone, I,J , Var,Err,Para,d1,d2,fileDEST)
% %  % ************************************  %  %
% %  % *********  By Yu-Hung Chang  *******  %  %
% %  % ***********  2016-Aug-29  **********  %  %
% %  % ************************************  %  %

IJ = I*J;
% %  % **** Err: lines between zones ******  %  %
% %  % **** Zone: Zones in one CFD case ***  %  %
% %  % **** Var: number of variables  *****  %  %
% %  % **** countNUM: number of snapshots *  %  %
% %  % **** CASE: number of Cases *********  %  %

 
fileID = fopen(fileDEST);
A = fread(fileID,'*single')';
fclose(fileID);

B = flipud(A'); clear A;
num = (IJ*Var+Err)*Zone;
B(num+1:end,:) =[];

DATA = reshape(B,(IJ*Var+Err),[]); clear B;
DATA(IJ*Var+1:end,:) =[];

OutPut = reshape(permute(reshape(DATA',Zone,IJ,Var),[2 1 3]),[],Var);

clear DATA
OutPut = flipud(OutPut);
OutPut = fliplr(OutPut);

OutPut((OutPut(:,1)> d1)|(OutPut(:,2)> d2)) = NaN;
OutPut2 =permute(reshape(OutPut,[IJ,Zone,Var]),[1,3,2]);clear OutPut;
idxToRemove = any(any(isnan(OutPut2),2),1);
OutPut2(:,:,idxToRemove) = [];
OutPut3 = reshape(permute(OutPut2,[2 1 3]),Var,[])'; clear OutPut2;

Results(:,1) = OutPut3(:,Para);
XYZ = OutPut3(:,1:3); 
end
