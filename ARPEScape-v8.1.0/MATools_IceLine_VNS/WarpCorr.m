function ACorr=WarpCorr(Angle,Energy,HV,ep,warpFitFile)
% ACorr=WarpCorr(Angle,Energy [,HV] [,ep] [,warpFitFile]) returns matrix ACorr of warping corrected angles 
% for the vector of angles Angle, vector of (kinetic or binding) energies Energy and pass energy ep. 
% If HV is omitted or empty, Energy is in kinetic energies; if a scalar or raw vector of photon energies, 
% Energy is in binding energies. If ep is omitted or empty, it is set to the default value based on the 
% input energy range. If the file of the warping fit coefficients warpFitFile is omitted or empty, it is 
% set to the default file WarpFit.dat.
% Note1: WarpCorr assumes the raw size of the data array as returned by ReadARPES without trimming
% Note2: WarpCorr applies to the data acquired in the snapshot mode but not in the scanning mode.
% Note3: WarpCorr keeps between the calls the correction matrices if the Angle, Energy and ep
% inputs are the same

disp('- Warping correction')

persistent iniChk WarpFit A1 A2 A3 A4 A5 A6 A7 A8 A9

% parameters
ePhi=4.5;
if nargin<3; HV=[]; end
if nargin<4; ep=[]; end
if nargin<5; warpFitFile='WarpFit.dat'; end;
% mode WAM/MAM/LAD/MAD = 1/2/3/4
aRange=Angle(end)-Angle(1);
mode=4; if aRange>12; mode=3; end; if aRange>18; mode=2; end; if aRange>24; mode=1; end;
% transformation to kinetic energy scale
if isempty(HV)
% - check against binding energy input
   if min(Energy)<4.5; disp('Warning: Binding energy scale requires photon energy input'); end
% - kinetic energy input    
   KEnergy=Energy;
else
   if min(Energy)<4.5  
% - binding energy input
      KEnergy=repmat(HV,length(Energy),1)+repmat(Energy,1,length(HV))-ePhi;
   else
% - kinetic energy input
      KEnergy=Energy;
   end
end
% default pass energy
eRange=Energy(end)-Energy(1);
if isempty(ep); ep=eRange/(2*0.066); end
% initialization if the first call or input energy or retarding ratio range changed
if isempty(iniChk)||norm(iniChk-[aRange eRange ep])>1e-3
% - load warping coefficients
   WarpFit=load(warpFitFile);
   ind=find(WarpFit(:,1)==mode); 
   if isempty(ind); disp('- Error: Angular mode not in WarpFit.dat'); ACorr=[]; return
   else WarpFit=WarpFit(ind,2:end);
   end
% - polynomial matrices
   [AM,EScM]=meshgrid(Angle,(Energy-mean(Energy))/ep);
   A1=AM; A2=AM.^2; A3=AM.^3; A4=AM.*EScM; A5=AM.^2.*EScM; 
   A6=AM.^3.*EScM; A7=AM.*EScM.^2; A8=AM.^2.*EScM.^2; A9=AM.^3.*EScM.^2;
% - initialization check
   iniChk=[aRange eRange ep];
end

% array of corrected angles
ACorr=zeros(size(KEnergy,1),length(Angle),size(KEnergy,2));
% - cycle over hv frames
for iS=1:size(KEnergy,2)
% - rr
rr=mean(KEnergy(:,iS))/ep;
% - interpolating the coefficients
Fit=interp1(WarpFit(:,1),WarpFit(:,2:end),rr,'pchip','extrap');
% - correction
ACorr(:,:,iS)=Fit(1)*A1+Fit(2)*A2+Fit(3)*A3+Fit(4)*A4+Fit(5)*A5+...
              Fit(6)*A6+Fit(7)*A7+Fit(8)*A8+Fit(9)*A9;
% - check corrected image
% ImData(ACorr,KEnergy,Data,'interp'); set(gca,'YDir','normal','TickDir','Out');
end