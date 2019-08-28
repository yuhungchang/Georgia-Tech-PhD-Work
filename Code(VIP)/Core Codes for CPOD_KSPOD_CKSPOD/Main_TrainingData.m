%% *****************************************
% ***************  2018-11-13 ***************
% ************** Yu-Hung Chang **************
% *******************************************
%% *****************************************

%% ===================== Loading raw data ===================== %%
Kt = Info.Kt; I = Info.I; J = Info.J;  Var = Info.Var;
Err = Info.Err; SKIP = Info.SKIP; Para = Info.Para;
d1 = Info.d1; d2 = Info.d2; List = Info.List;
Ln = Info.Ln;  Rn = Info.Rn; CASE = Info.CASE; SNAP = Info.SNAP;
DATA = cell(CASE,1);  XYZ0 = cell(CASE,1);
clc
fprintf('\n *********** Loading raw data *********** \n');
tic

if Core <=CASE
    CC = Core;     
else
    CC = CASE;     
end

parfor (Lp = 1:CASE,CC)  %Case number
    tic
    k1 = int2str(List(Lp));
    Folder = ['../case',k1,'_plt'];
%     Folder  = fullfile(data2/EmulationPLT,DES);
    [DATA{Lp},XYZ0{Lp}] = LoadingPLT(Kt, I,J , Var,Err,SKIP,Para,d1,d2,Folder);
    DATA{Lp} = double(DATA{Lp});
    XYZ0{Lp} = double(XYZ0{Lp});
    fprintf('\n *********** Complete loading Case %s for %d files *********** \n',k1,SNAP);
    toc
end
End = toc;  Min = floor(End/60); Sec = End-Min*60;
% delete(gcp);
XYZ = XYZ0{1}; clear XYZ0
XYZ = double(XYZ);  clc;
fprintf('\n --------- Total time for loading %d PLT files: %d min %5.2f sec.  -------------------- \n',CASE*SNAP,Min,Sec);
fprintf('\n --------- End up loading data : %d Cases (%d files per case).  -------------------- \n',CASE,SNAP);
 