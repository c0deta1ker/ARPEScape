InputFile='X:\CommiData\2015\08_2015\CurvatureWAM\WAMC_Ep100.h5'; mode=1;  % curvature test
% InputFile='X:\CommiData\2013\07_2013\WAM_Uc_Ep=100.sp2'; mode=1;  % curvature test
% InputFile='X:\CommiData\2014\02_2014\MAM_Uc_Ep=60.sp2'; mode=2; % curvature test
% InputFile='X:\CommiData\2014\01_2014\WarpMAM_Ek600_Ep200.sp2'; mode=2; % warping test

% SL2 version
% [~,Energy,~,~,ep,~]=ReadARPES(InputFile); ek=mean(Energy); 
% Data=ReaderSP2(InputFile); Data=(Data.raw)'; Data=Data(1:2:end,1:2:end);

% Prodigy version
[~,Energy,~,Data,ep]=ReadARPES([IniDir Files{iFile}]); ek=mean(Energy);

TransformFile='Transform_Raw.dat';
[ARange,ERange,DataT]=Transform_Raw(Data,ek,mode,ep,TransformFile);

ImData(linspace(ARange(1),ARange(end),size(DataT,2)), linspace(ERange(1),ERange(2),size(DataT,1)), DataT,'flat')