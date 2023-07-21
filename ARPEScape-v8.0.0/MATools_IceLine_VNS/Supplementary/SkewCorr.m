function DataCorr=SkewCorr(Angle,Energy,Data,AngleC,EnergyC,DataC,eLoc,feat)
% SkewCorr.m removes skewing of Data in energy scale as function of angle, which may result from 
% certain irreproducibility of the analyzer slit setting. The algorithm uses 2-order polynomial fit 
% of the Fermi edge or spectral peak in the calibration data, averaged within 20 angular bins, and 
% AlignEF.m with default energy window and pre-smoothing parameters.
% Inputs: Angle (1D-row), Energy (1D-column) and Data (2D or 3D) are the input data arrays; 
% AngleC (1D-row), EnergyC (1D-column) and DataC (2D) are the calibration data arrays; optional eLoc
% (default eLoc=0) specifies approximate energy of the calibration spectral structure, and optional 
% feat = 'EF' or 'Peak' (case-insensitive, default 'EF') the type of this structure. 
% Outputs: DataCorr as a function of the same Angle and Energy

% calibration data and parameters
%[AngleC,EnergyC,~,DataC,~,~]=ReadARPES(CalibFile); % the previous version 
nA=size(DataC,2); nBin=20;
if nargin<7; eLoc=[]; end; if isempty(eLoc); eLoc=0; end
if nargin<8; feat=[]; end; if isempty(feat); feat='EF'; end;

% alignment and fit of the calibration data
% - global EF alignment to help reliability of the local alignment
[EnergyC,~,Fail]=AlignEF(DataC,EnergyC,eLoc,[],[],feat);
if ~isempty(Fail); DataCorr=[]; disp('Global EF alignment error'); return; end
% - local EF alignment
width=floor(nA/nBin);
n1=round((nA-nBin*width)/2);
A=[]; EF=[];
for j=n1:width:n1+width*(nBin-1);
   [~,ef,Fail]=AlignEF(DataC(:,j:j+width),EnergyC,0,[],[],feat);
   if isempty(Fail); 
      A=[A (AngleC(floor(j+width/2))+AngleC(ceil(j+width/2)))/2]; EF=[EF ef]; 
   end
end
% - polynomial fit 
P=polyfit(A,EF,2);
% figure; plot(A,EF,'b',AngleC,polyval(P,AngleC),'r')
% - removal of constant offset
P(end)=0;

% input data alignment
% - interpolation matrices
AngleM=repmat(Angle,size(Data,1),1);
EnergyCorr=repmat(Energy,1,size(Data,2))+repmat(polyval(P,Angle),size(Data,1),1);
% - cycle through scans 
DataCorr=Data; for iScan=1:size(Data,3)
% - interpolating the corrected data
   DataCorr(:,:,iScan)=interp2(Angle,Energy,Data(:,:,iScan),AngleM,EnergyCorr,'pchip');
end