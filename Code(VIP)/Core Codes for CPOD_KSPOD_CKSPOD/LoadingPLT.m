%% *****************************************
% ***************  2018-11-13 ***************
% ************** Yu-Hung Chang **************
% *******************************************
%% *****************************************
function [Results,XYZ] = LoadingPLT(Zone, I,J , Var,Err,SKIP,Para,d1,d2,Folder)
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
% Snap = countNUM-countStart+1;
LL = (1:(SKIP+1):1000)'; 
% Results = zeros(IJ*Zone,length(LL));
Loop =1; point =1;


for PP= 1:length(LL)
    H2 = LL(PP);
    k2 = int2str(LL(PP));
    if H2<=9
        filename=['plot000',k2,'.plt'];
    elseif (H2>=10 && H2<=99)
        filename=['plot00',k2,'.plt'];
    elseif (H2>=100 && H2<=999)
        filename=['plot0',k2,'.plt'];
    else
        filename=['plot',k2,'.plt'];
    end
    
    fileDEST  = fullfile(Folder,filename);
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
    
    Results(:,point) = OutPut3(:,Para);
    point = point+1;
    if Loop ==1
        XYZ = OutPut3(:,1:3);
        
        Loop= Loop+1;
    end
    clear OutPut
    %         clearvars -except Err Var KK IJ OLD Zone Results H1 H2 CASE countNUM
end

end
