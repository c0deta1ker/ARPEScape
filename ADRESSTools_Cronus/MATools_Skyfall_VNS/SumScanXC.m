function [DataXC,OffsE]=SumScanXC(Energy,Data,maxLagE,Win)
% [DataXC,OffsE] = SumScanXC(Energy,Data,maxLagE[,WinE]) summs up the separate scan frames from 3D 
% array Data with energy offsets calculated by cross-correlation within the lag interval % +-maxLagE. 
% The cross-correlation energy window is Win (default 80% of Energy) and angular window is fixed to 
% default 90%. The returned offsets OffsE are discrete on energy grid without interpolation.
% Note: The cross-correlation algorithm subtracts linear background fixed at the end points of the 
% energy window. For better accuracy choose the window ends at valleys between the spectral peaks 
% and not far above EF.  
% Ver. 17.06.2017
disp('- Correlating scans')
% check input
if size(Data,3)<2; DataXC=Data; OffsE=0; return; end
% parameters
nE=length(Energy); dE=mean(diff(Energy)); nA=size(Data,2); nS=size(Data,3);
maxLag=round(maxLagE/dE); Lag=-1*maxLag:maxLag;
% XC-window
% - default window of 80% in energy and 90% in angle
EWin=(1+round(0.1*nE)):round(0.9*nE);
AWin=(1+round(0.05*nA)):round(0.95*nA); nAWin=length(AWin);
% - energy window
if nargin>3;
   Win=sort(Win);
   EWin=round([max(1,interp1(Energy,1:nE,Win(1))):min(nE,interp1(Energy,1:nE,Win(2)))]);
end
nEWin=length(EWin);

% data window
DWin=Data(EWin(1):EWin(end),AWin(1):AWin(end),:); 
X=repmat(1:nAWin,nEWin,1); Y=repmat((1:nEWin)',1,nAWin);
% figure; surf(X,Y,DWin(:,:,1),'edgecolor','none');
% suppress the window edge effects by 1D linear fit of the edge values
% - window edge values averaged over 3 end points
D1=sum(DWin(1:3,:,:),1)/3; D2=sum(DWin(end-2:end,:,:),1)/3;
for iS=1:nS
  D1S=D1(1,:,iS); D2S=D2(1,:,iS);
  Bkg=repmat(D1S,nEWin,1)+repmat((D2S-D1S)./(Y(end,:)-Y(1,:)),nEWin,1).*Y;
%  figure; surf(X,Y,Bkg,'edgecolor','none');
  DWin(:,:,iS)=DWin(:,:,iS)-Bkg; 
%  figure; surf(X,Y,DWin(:,:,iS),'edgecolor','none');
end

% padded reference
Ref=DWin(:,:,1); RPad=[zeros(maxLag,nAWin);Ref;zeros(maxLag,nAWin)];
% cycle over scans
OffsE=0; DataXC=Data(:,:,1)/nS;
for iS=2:nS
% - cycle over lags
   XC=[]; for lag=Lag
      DPad=[zeros(maxLag+lag,nAWin);DWin(:,:,iS);zeros(maxLag-lag,nAWin)];
      xc=sum(sum(RPad.*DPad));
%/(size(Ref,1)-abs(lag)); % Note: Unbiased correlation may fail for flat-topped functions like y=x^2 (see TestXC.m)
      XC=[XC xc];
   end
% figure; plot(dE*Lag,XC);
% - offset at maximal correlation 
   [~,indMax]=max(XC); offs=Lag(indMax); offsE=offs*dE; OffsE=[OffsE; offsE];
   disp(['Scan = ' num2str(iS) ': Offset = ' num2str(offsE)])
% - appending/truncation of the correlated images
%   Note: negative arguments of ones(n,m) return empty array 
Data(:,:,iS)=[NaN*ones(offs,nA); ...
                 Data(max(1,1-offs):min(nE,nE-offs),:,iS); ...
                 NaN*ones(-1*offs,nA)];
% - summation normalized by nS
   DataXC=DataXC+Data(:,:,iS)/nS;
end