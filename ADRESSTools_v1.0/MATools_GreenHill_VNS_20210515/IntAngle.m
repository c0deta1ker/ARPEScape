function [DataInt,EInt]=IntAngle(Data,ACorr,ECorr,Interval)
% [DataInt,EInt]=IntAngle(Data,ACorr,ECorr [,Interval]) integrates 2D/3D
% array Data over 1D/2D array ACorr of angles within Interval (default 
% is 90% of the ACorr range). The outputs are respectively 1D/2D integrated 
% data array DataInt and 1D energy vector EInt
% Rev. 07.07.2019 

disp('- Angle integrating')
% parameters
% - angular interval 
if nargin<4; Interval=0.9*[max(ACorr(:,1)) min(ACorr(:,end))]; end
% - dimensions
nE=size(Data,1); nA=size(Data,2); nS=size(Data,3);
% - extend energy array to the data array
if size(ECorr,2)==1; ECorr=repmat(ECorr,1,nA); end
if nS>1&&size(ECorr,3)==1; ECorr=repmat(ECorr,1,1,nS); end
% - reference energy grid
if nA==1&&nS==1; EInt=ECorr;
else eMin=max(max(ECorr(1,:,:))); eMax=min(min(ECorr(end,:,:))); EInt=linspace(eMin,eMax,nE)'; 
end
% data within the angular interval/normalization
if size(ACorr,1)==1
% - 1D-ACorr    
   Data=Data(:,Interval(1)<=ACorr&ACorr<=Interval(2),:); 
   nA=size(Data,2); Data=Data/nA;
else
% - 2D-ACorr
   Out=find(Interval(1)>ACorr|ACorr>Interval(2)); Data(Out)=0;
   Supp=ones(size(Data)); Supp(Out)=0;
   Weight=sum(Supp,2); Weight=repmat(Weight,1,nA);
   Data=Data.*Supp./Weight;
end

% cycle over scans
DataInt=zeros(nE,nS);
for iS=1:nS
% cycle over angles in triads for speed
   IntScan=zeros(size(EInt));
   for iA=2:3:nA-1
% interpolation on the reference grid
      if nS>1
         DCorr=interp1(ECorr(:,iA,iS),sum(Data(:,iA-1:iA+1,iS),2),EInt);
      else
         if nA>1
            DCorr=interp1(ECorr(:,iA),sum(Data(:,iA-1:iA+1),2),EInt);
         else DCorr=interp1(ECorr,sum(Data(:,iA-1:iA+1),2),EInt);
         end
      end
% angle integration
      IntScan=IntScan+DCorr;
   end
% check
%   figure; plot(EInt,DD,EInt,DDS);
% output and normalization
   DataInt(:,iS)=IntScan;
end