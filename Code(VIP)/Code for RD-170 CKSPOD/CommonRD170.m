function  [Common,commonXYZ,Label,Node] = CommonRD170(CASE,SNAP,xyzU,xyzD,DataU,DataD,MeshFile,FileDir,delimiterIn,headerlinesIn)
% % % ************************************  % % %
% % % *********  By Yu-Hung Chang  *******  % % %
% % % ***********  2017-March-26  ********  % % %
% % % ************************************  % % %
tic
for H=1:3 
    filename = char(MeshFile(H));
    fileDEST  = fullfile(FileDir,filename);    
    
    mesh0=importdata(fileDEST,delimiterIn,headerlinesIn);
%     mesh0=importdata('D:\NewEmulationData\RD-170Group\Mesh_zone57-96_merge.dat',delimiterIn,headerlinesIn);
    mesh1 = struct2cell(mesh0);
    mesh2 = cell2mat(mesh1(1)); 
    Mesh{H} = mesh2(~any(isnan(mesh2),2),1:2);
%     eval(['Mesh' num2str(H) ' = mesh2(~any(isnan(mesh2),2),:);']);
    
    mesh3 = mesh2(any(isnan(mesh2),2),:);
    Node{H}  = mesh3(:,~any(isnan(mesh2),1)); 
    clear mesh0 mesh1 mesh2 filename fileDEST mesh3
end
t1= toc;
fprintf('\n *********** Complete loading common grid mesh, time = %5.2f sec *********** \n',t1);
  

%% Labeling
Mesh = Mesh';
L1 = ones(size(Mesh{1},1),1);
L2 = ones(size(Mesh{2},1),1)*2;
L3 = ones(size(Mesh{3},1),1)*3;
% L4 = ones(size(Mesh{4},1),1)*4;

Label = [L1;L3]; 
commonXYZA = Mesh{1}; %==== Upstream Only =
commonXYZB = Mesh{3}; %==== Downstream Only =
commonXYZ = [commonXYZA;commonXYZB]; 

% warning('off','all')
% warning

%% Common Grid Interpolation
Common = cell(CASE,1);
XYZ = cell(CASE,1);
Data = cell(CASE,1);

Time = zeros(CASE,1);
tic
parfor (Lp=1:CASE,6)
    tic     
    XYZ{Lp} = [xyzU{Lp};xyzD{Lp}];
    Data{Lp} = [DataU{Lp};DataD{Lp}];   
    
%     XYZ{Lp} = [xyzD{Lp}];
%     Data{Lp} = [DataD{Lp}]; 
    
    for fil =1:SNAP
        F = scatteredInterpolant(XYZ{Lp}(:,1),XYZ{Lp}(:,2),Data{Lp}(:,fil));
        Area = F(commonXYZ(:,1),commonXYZ(:,2));
        Common{Lp}(:,fil) = Area;         
    end
    Time(Lp) = toc;
    fprintf('\n *********** Complete common grid for case %d, time = %5.2f sec *********** \n',Lp,Time(Lp));
end
clc;  
clear XYZ Data
  
totalT = sum(Time)+t1;
fprintf('\n *********** Complete common grid process for %d cases are done! Total time = %5.2f sec *********** \n',CASE,totalT);


   
end


