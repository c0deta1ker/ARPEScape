% Fitting the warping pattern produced by slit array. The input files are
% either sp2 (kinetic energy scale) or hdf5 (kinetic or binding energy 
% scale). The warping is a function of angle and energy (the latter 
% normalized to the [0 1] detector length)

% parameters
check=1;
ePhi=4.5;
CurveCorrFile='CurveFit_Raw.dat';
aTrim=0.5;

% warping calibration 15.10.2018 after the reparation at SPECS/ setting Icoil=-2mA
% % angular offset determination
% Prefix='X:\CommiData\2018\15.10.2018\MAMW_Ek600_CentreLine_Ep'; FileN=100; Ext='.h5';
% mode=2; nFit=1; ECut=(linspace(-0.85,0.85,20))'; XOrig=0;
% % - WAM
% Prefix='X:\CommiData\2018\25.09.2018\WAMW_Ek350_Ep'; Ext='.h5';
% FileN=[30:5:100 110:10:230]; %FileN=30;
% mode=1; nFit=21; ECut=(linspace(-0.875,0.875,20))'; XOrig=0.028;
% % - MAM
% Prefix='X:\CommiData\2018\15.10.2018\MAMW_Ek600_Ep'; Ext='.h5';
% FileN=[45:5:100 110:10:210 230:20:350];  % FileN=100;
% mode=2; nFit=17; ECut=(linspace(-0.85,0.85,20))'; XOrig=0;
% % - LAD
% Prefix='X:\CommiData\2018\25.09.2018\LADW_Ek800_Ep'; Ext='.h5';
% FileN=[35:5:100 110:10:210 230:20:350];  % FileN=35:5:90;
% mode=3; nFit=11; ECut=(linspace(-0.925,0.825,20))'; XOrig=0;
% % Note: Discontinuity in determined peaks breaks the fit

% Curvature correction after the upgrade to PHOIBOS225 in 2022
% - MAM
Prefix='X:\ADRESS_e18633\CommiData\2022\12-14.08\WarpingCalibration\MAMW_Ek600_Ep'; Ext='.h5';
FileN=[30:5:100 110:10:210 230:20:350];  % FileN=100;
mode=2; nFit=19; nShow=21; ECut=(linspace(-0.9,0.9,20))'; XOrig=0;

% energy correction table
fid=fopen(CurveCorrFile);
TabE=[]; while ~feof(fid)
   DLine=fgetl(fid); DLine=sscanf(DLine,'%f'); DLine=DLine';
   if ~isempty(DLine)
      if DLine(1)==100+mode; TabE=[TabE; DLine(2:end)]; end
   end   
end
fclose(fid);

% file loop
WarpFit=[]; for iFile=1:size(FileN,2)
File=[Prefix num2str(FileN(iFile)) Ext];
[~,Energy,~,Data,Ep]=ReadARPES(File);
% disp(['File = ' File])
% % - angular ROI and offset
% nA=size(Data,2); xOffset=nA/125; Data=Data(:,round(nA/8+xOffset):round(7*nA/8+xOffset));
% - parameters
nA=size(Data,2); nE=size(Data,1);
% - coordinates normalized to the [-1 1] interval
ASc=linspace(-1,1,nA); ESc=linspace(-1,1,nE); [AScM,EScM]=meshgrid(ASc,ESc);
% % - hv and kinetic energies from hdf5 file
% try pos1=strfind(Note,'hv'); pos2=strfind(Note,'Pol'); hv=str2num(Note(pos1+10:pos2-1)); end
% if mean(Energy)<0; Energy=Energy+hv-ePhi; end;
% - RR
rr=mean(Energy)/Ep;
% - check input
if check==1; figure;
   subplot(3,1,1); ImData(ASc(1:2:end),ESc(1:2:end),Data(1:2:end,1:2:end),'interp');
   title(File,'Interpreter','none','fontweight','normal');
end

% curvature correction
FitE=interp1(TabE(:,1),TabE(:,2:end),Ep,'pchip','extrap');
% - NOTE: If Transform_Raw attempts to correct the curvature and warping simultaneously in the 
% same interp2, the following line shall be commented out. However, in some way the simultaneous 
% correction distorts the curvature by some kind of crosstalk. Therefore, this code assumes that 
% the curvature is corrected in one interp2 and the warping in another interp2. 
DE=polyval([FitE 0],ASc); ECorr=EScM+repmat(DE,nE,1); 
Data=interp2(AScM,EScM,Data,AScM,ECorr,'cubic'); Data(isnan(Data))=0;
% - check
% figure; ImData(ASc(1:2:end,1:2:end),ESc(1:2:end,1:2:end),Data(1:2:end,1:2:end),'interp');

% loop through the cuts
XMaxM=[]; XMaxC=[]; AngM=[];
for iCut=1:length(ECut)
% - selecting the slice and summing up    
   dCut=mean(diff(ECut));
   A=Data(ESc>ECut(iCut)-dCut/2&ESc<ECut(iCut)+dCut/2,:); A=sum(A,1)/size(A,1);
% - trim below half-height to remove the stray maxima
   A=A-aTrim*mean(A); A(A<0)=0;
% - Gaussian smoothing; 'extrap' to suppress stray edge maxima    
   hwPts=0.0075*(nA-1); A=VGaco(A,hwPts,3*hwPts,'extrap');
% - maxima   
   [~,XMax]=ExLoc(ASc,A);
% - central maximum
   [~,indC]=min(abs(XMax-XOrig)); XMaxC=[XMaxC;XMax(indC)];
% - number of maxima below and above the central one
% nXMaxLo=indC-length(XMax(XMax<XMax(indC))); nXMaxHi=indC+length(XMax(XMax>XMax(indC)));
nXMaxLo=length(XMax(XMax<XMax(indC))); nXMaxHi=length(XMax(XMax>XMax(indC)));
% - trimming the number of maxima
% nXMaxLo=max(nXMaxLo,1); nXMaxHi=min(nXMaxHi,nFit);
nXMaxLo=min(nXMaxLo,(nFit-1)/2); nXMaxHi=min(nXMaxHi,(nFit-1)/2);
% XMax=XMax(nXMaxLo:nXMaxHi);
% XMax=XMax(indC-(nFit-1)/2:indC+(nFit-1)/2);
XMax=XMax(indC-nXMaxLo:indC+nXMaxHi);
indC=1+nXMaxLo;
% - mapping the maxima onto 1:nFit, with the missing maxima set to NaN
nOffset=1+(nFit-1)/2-indC;
XMaxApp=NaN*ones(1,nFit); XMaxApp(1+nOffset:nOffset+length(XMax))=XMax; XMax=XMaxApp;   
% - append to matrix 
   XMaxM=[XMaxM;XMax];
% - angular grid
   dA=1; AngM=[AngM;dA*((1:length(XMax)))-(nFit-1)/2-1];
end
% check
if check==1; subplot(3,1,2); plot(XMaxM,ECut,'-'); end; ylim([-1 1])

% alignment of the central line
% - fitting
LineFit=polyfit(ECut,XMaxC,1); XMaxFit=polyval(LineFit,ECut);
% - alignment
XMaxM=XMaxM-repmat(XMaxFit,1,size(XMaxM,2));
% - check
% subplot(3,1,3); plot(XMaxM,ECut,'-'); ylim([-1 1])

% linear fit
% - flatten and offset the arrays
XMaxF=XMaxM(:);
EM=repmat(ECut,1,size(XMaxM,2)); EF=EM(:);
% - NOTE: If Transform_Raw attempts to correct the curvature and warping simultaneously
% in the same interp2 (see the above % curvature correction section), the following line 
% shall be introduced here: DE=polyval([FitE 0],XMaxF); EF=EF-DE;
AngF=AngM(:);
% - remove NaNs
IRem=isnan(XMaxF); XMaxF(IRem)=[]; EF(IRem)=[]; AngF(IRem)=[];
% - fit (of the order 4 in angle and 3 in energy as optimized with actual data)
FString='X X.^2 X.^3 X.^4 X.*Y X.^2.*Y X.^3.*Y X.^4.*Y X.*Y.^2 X.^2.*Y.^2 X.^3.*Y.^2 X.^4.*Y.^2';
% FString='X X.^2 X.^3 X.^4 X.^5 X.*Y X.^2.*Y X.^3.*Y X.^4.*Y X.^5.*Y X.*Y.^2 X.^2.*Y.^2 X.^3.*Y.^2 X.^4.*Y.^2 X.^5.*Y.^2';
Fit=LinFit(AngF,EF,XMaxF,FString);
% - check the corrected image:
if check==1
% - - coordinate grids
   AGrid=linspace(-1*(nShow-1)/2-dA/2,(nShow-1)/2+dA/2,nA); [AGridM,~]=meshgrid(AGrid,ESc);
% - - back transformation of physical angle coordinates to pixel coordinates
   ACorr=Fit(1)*AGridM+Fit(2)*AGridM.^2+Fit(3)*AGridM.^3+Fit(4)*AGridM.^4+...
         Fit(5)*AGridM.*EScM+Fit(6)*AGridM.^2.*EScM+Fit(7)*AGridM.^3.*EScM+Fit(8)*AGridM.^4.*EScM+...
         Fit(9)*AGridM.*EScM.^2+Fit(10)*AGridM.^2.*EScM.^2+Fit(11)*AGridM.^3.*EScM.^2+Fit(12)*AGridM.^4.*EScM.^2;
%    ACorr=Fit(1)*AGridM+Fit(2)*AGridM.^2+Fit(3)*AGridM.^3+Fit(4)*AGridM.^4+Fit(5)*AGridM.^5+...
%          Fit(6)*AGridM.*EScM+Fit(7)*AGridM.^2.*EScM+Fit(8)*AGridM.^3.*EScM+Fit(9)*AGridM.^4.*EScM+Fit(10)*AGridM.^5.*EScM+...
%          Fit(11)*AGridM.*EScM.^2+Fit(12)*AGridM.^2.*EScM.^2+Fit(13)*AGridM.^3.*EScM.^2+Fit(14)*AGridM.^4.*EScM.^2+Fit(15)*AGridM.^5.*EScM.^2;   
% - - data array interpolated on the back transformed pixel coordinates
   DataCorr=interp2(AScM,EScM,Data,ACorr,EScM,'cubic');
% - check
   subplot(3,1,3); ImData(AGrid(1:2:end,1:2:end),ESc(1:2:end,1:2:end),DataCorr(1:2:end,1:2:end),'interp');
end
% - fit array (mode offset by 10)
WarpFit=[WarpFit;[mode+10 rr Fit]];
% end of the file loop
end

% save the fit data
save WarpFit_Raw.dat WarpFit -ascii