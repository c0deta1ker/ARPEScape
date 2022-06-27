function DataCorr=PolyAlignEF(Angle,Energy,Data,CalibFile)
% PolyCorrEF.m aligns Data in energy scale relative to EF varying as function of angle. The 
% algorithm uses 3-order polynomial fit of the Fermi edge, averaged within 20 angular bins.
% Inputs: Angle (1D-row), Energy (1D-column) and Data (2D or 3D) are the input array; CalibFile
% is the EF calibration file. Outputs: DataCorr as a function of the same Angle and Energy

% calibration data and parameters
[AngleC,EnergyC,~,DataC,~,~]=ReadARPES(CalibFile);
nA=size(DataC,2); nBin=20;
% - global EF alignment to help reliability of the local alignment
[EnergyC,~,Fail]=AlignEF(DataC,EnergyC,0);
if ~isempty(Fail); DataCorr=[]; disp('Global EF alignment error'); return; end
% - local EF alignment
width=floor(nA/nBin);
n1=round((nA-nBin*width)/2);
A=[]; EF=[];
for j=n1:width:n1+width*(nBin-1);
   [~,ef,Fail]=AlignEF(DataC(:,j:j+width),EnergyC,0);
   if isempty(Fail); 
      A=[A (AngleC(floor(j+width/2))+AngleC(ceil(j+width/2)))/2]; EF=[EF ef]; 
   end
end
% - polynomial fit 
P=polyfit(A,EF,3);
% figure; plot(A,EF,'b',AngleC,polyval(P,AngleC),'r')

% input data alignment
% - interpolation matrices
AngleM=repmat(Angle,size(Data,1),1);
EnergyCorr=repmat(Energy,1,size(Data,2))+repmat(polyval(P,Angle),size(Data,1),1);
% - cycle through scans 
DataCorr=Data; for iScan=1:size(Data,3)
% - interpolating the corrected data
   DataCorr(:,:,iScan)=interp2(Angle,Energy,Data(:,:,iScan),AngleM,EnergyCorr,'pchip');
end