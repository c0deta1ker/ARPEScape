% Fitting the curvature along the detector

% Curvature before the slit upgrade in September 2015
% % WAM, Ag3d, Uconv=-0.1*Ep
% IniDir='X:\CommiData\2013\07_2013\';
% Files={'WAM_Uc_Ep=30.sp2' 'WAM_Uc_Ep=35.sp2' 'WAM_Uc_Ep=40.sp2' 'WAM_Uc_Ep=50.sp2' 'WAM_Uc_Ep=60.sp2' ...
%        'WAM_Uc_Ep=70.sp2' 'WAM_Uc_Ep=80.sp2' 'WAM_Uc_Ep=90.sp2' 'WAM_Uc_Ep=100.sp2' 'WAM_Uc_Ep=110.sp2' ...
%        'WAM_Uc_Ep=120.sp2' 'WAM_Uc_Ep=130.sp2' 'WAM_Uc_Ep=150.sp2' 'WAM_Uc_Ep=170.sp2' 'WAM_Uc_Ep=190.sp2' ...
%        'WAM_Uc_Ep=210.sp2' 'WAM_Uc_Ep=230.sp2' };
% mode=1; aLo=-0.9; aHi=0.9; nPeak=1;
% % MAM, Au 4f, Uconv=-0.1*Ep
% IniDir='X:\CommiData\2014\02_2014\';
% Files={'MAM_Uc_Ep=50.sp2' 'MAM_Uc_Ep=60.sp2' 'MAM_Uc_Ep=70.sp2' 'MAM_Uc_Ep=80.sp2' ...
%        'MAM_Uc_Ep=90.sp2' 'MAM_Uc_Ep=100.sp2' 'MAM_Uc_Ep=110.sp2' 'MAM_Uc_Ep=120.sp2' 'MAM_Uc_Ep=130.sp2' ...
%        'MAM_Uc_Ep=150.sp2' 'MAM_Uc_Ep=170.sp2' 'MAM_Uc_Ep=190.sp2' 'MAM_Uc_Ep=210.sp2' 'MAM_Uc_Ep=240.sp2' ...
%        'MAM_Uc_Ep=280.sp2' 'MAM_Uc_Ep=320.sp2' 'MAM_Uc_Ep=360.sp2' 'MAM_Uc_Ep=400.sp2'};
% % MAM_Uc_Ep=55.sp2 yields incorrect peak position     
% mode=2; aLo=-0.9; aHi=0.875; nPeak=1;
% % LAD, Au 4f, Uconv=-0.1*Ep
% IniDir='X:\CommiData\2014\06_2014\';
% Files={'LAD_Uc_Ep20.sp2' 'LAD_Uc_Ep25.sp2' 'LAD_Uc_Ep30.sp2' 'LAD_Uc_Ep35.sp2' 'LAD_Uc_Ep40.sp2' ...
%        'LAD_Uc_Ep50.sp2' 'LAD_Uc_Ep60.sp2' 'LAD_Uc_Ep80.sp2' 'LAD_Uc_Ep100.sp2' 'LAD_Uc_Ep120.sp2' 'LAD_Uc_Ep140.sp2' ...
%        'LAD_Uc_Ep170.sp2' 'LAD_Uc_Ep200.sp2' 'LAD_Uc_Ep250.sp2' 'LAD_Uc_Ep300.sp2' 'LAD_Uc_Ep350.sp2' 'LAD_Uc_Ep380.sp2'};
% % LAD_Uc_Ep16.sp2 is unreliable   
% % mode=3; aLo=-0.825; aHi=0.8; nPeak=1; % iFile<3
% mode=3; aLo=-0.9; aHi=0.875; nPeak=1; % iFile>3

% Curvature after the slit upgrade in September 2015
% IniDir='X:\CommiData\2015\05_10_2015\'; Suffix='.h5';
% % WAM, Au 4f, Ek=350 eV
% Prefix='WAMC_Ek350_Ep';
% FileN={'35' '40' '45' '50' '60' '70' '80' '90' '100' '110' '130' '150' '170' '190' '210' '230'};
% mode=1; aLo=-0.9; aHi=0.88; nPeak=1;
% % MAM, Au 4f, Ek=600 eV
% Prefix='MAMC_Ek600_Ep';
% FileN={'50' '55' '60' '65' '70' '75' '80' '90' '100' '110' '120' '130' '140' '150' ...
%           '170' '190' '210' '230' '270' '310' '350'};
% mode=2; aLo=-0.9; aHi=0.9; nPeak=1;
% % LAD, Au 4f, Ek=800 eV
% Prefix='LADC_Ek800_Ep';
% FileN={'35' '40' '45' '50' '55' '60' '65' '70' '75' '80' '90' '100' '110' '120' '130' '140' ...
%        '150' '170' '190' '210' '230' '270' '310' '350'};
% % '30' contains invalid data   
% mode=3; aLo=-0.88; aHi=0.85; nPeak=1;

% Curvature after the upgrade to PHOIBOS225 in 2022
% - MAM, Au 4f, Uconv=-0.1*Ep
IniDir='X:\ADRESS_e18633\CommiData\2022\12-14.08\CurvatureCalibration\';
Prefix='MAMC_Ep';
FileN={ '32_5' '35' '40' '45' '50' '55' '60' '65' '70' '75' '80' '90' '100' '110' '120' '130'... 
           '140' '150' '170' '190' '210'  '230' '270'  '310' '350' '400'};
mode=2; aLo=-1; aHi=1; nPeak=1;

% input
CurveFit=[]; for iFile=1:size(FileN,2)
File=[IniDir Prefix FileN{iFile} '.h5'];
[~,~,~,Data,ep]=ReadARPES(File);
%Data=ReaderSP2([IniDir Files{iFile}]); Data=(Data.raw)';

% % angular ROI used in Transform_Raw (comment out for PHOIBOS225)
% nA=size(Data,2); Data=Data(:,round(nA/8):round(7*nA/8));

% dimensions
nA=size(Data,2); Angle=linspace(-1,1,nA);
nE=size(Data,1); Energy=linspace(-1,1,nE);
% figure; ImData(Angle,Energy,Data)

% % mode WAM/MAM/LAD/MAD = 1/2/3/4
% aRange=Angle(end)-Angle(1);
% mode=4; if aRange>12; mode=3; end; if aRange>18; mode=2; end; if aRange>24; mode=1; end;

% angular ranger of the fit
ARange=find(Angle>=aLo&Angle<=aHi);
Angle=Angle(ARange); Data=Data(:,ARange);
% figure; ImData(Angle,Energy,Data)

% data reduction to triades
Angle=(Angle(1:3:end-2)+Angle(2:3:end-1)+Angle(3:3:end))/3;
Data=(Data(:,1:3:end-2)+Data(:,2:3:end-1)+Data(:,3:3:end))/3;

% evaluation of the extremes
EMax=[]; for iA=1:length(Angle)
   A=Data(:,iA);
% - smoothing
   hw=20; A=Gaco1(A,hw);
% - data trimming at half-height
   A=A-0.5*max(A); A(A<0)=0; 
% - extremes   
   [~,XMax]=ExLoc(Energy,A);
% - identification of the main extremes   
%   YMax=interp1(Energy,A,XMax); [~,Ind]=sort(YMax); XMax=XMax(Ind);   
   XMax=XMax(end);
% - nPeak largest peaks    
   XMax=XMax(end-(nPeak-1):end); eMax=mean(XMax);
   EMax=[EMax;eMax];
end

% fitting

% for the burn out MCP from PHOIBOS150   
EMaxTrim=EMax; EMaxTrim((Angle>-0.78&Angle<-0.67)|(Angle>0.66&Angle<0.82))=[];
AngleTrim=Angle; AngleTrim((Angle>-0.78&Angle<-0.67)|(Angle>0.66&Angle<0.82))=[];

% - fitting
% curveFit=polyfit(Angle',EMax,5);

curveFit=polyfit(AngleTrim',EMaxTrim,5); % for the burn out MCP from PHOIBOS150   

disp(num2str(curveFit));
% - running info
figure; plot(Angle,EMax,Angle,polyval(curveFit,Angle))
% output array (correction at zero angle set to zero)
CurveFit=[CurveFit; [mode+100 ep curveFit(1:end-1)]];
end

save CurveFit_Raw.dat CurveFit -ascii