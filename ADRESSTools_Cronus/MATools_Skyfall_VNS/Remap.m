function [XR,YR,DataR]=Remap(XM,YM,Data)
% [XR,YR,DataR]=Remap(XM,YM,Data) remaps the data array Data defined on the raw 
% vector/2D/3D array XM of x-coordinates (usually the corrected angles ACorr) and 
% column vector/2D/3D array YM of y-coordinates (usually aligned energies EAlign)
% onto DataR defined on the rectangular grid formed by the equidistant vectors XR and YR
% Ver. 8 Jan 2021

disp('- Re-mapping')

% check if XM and YM are already equidistant vectors
if ndims(XM)==1&&ndims(YM)==1
   if all([diff(diff(XM))==0 diff(diff(XM))]==0) 
      disp('- No remapping of the output data because the input coordinates are already equidistant vectors'); 
      XR=XM; YR=YM; DataR=Data; return 
   end
end

% parameters
nA=size(XM,2); nA3=size(XM,3);
nE=size(YM,1); nE3=size(YM,3); 
nS=size(Data,3);

% extending 1D arrays to 2D
if size(XM,1)==1; XM=repmat(XM,nE,1); end
if size(YM,2)==1; YM=repmat(YM,1,nA); end
% output grid
XR=linspace(min(min(min(XM))),max(max(max(XM))),nA);
YR=linspace(min(min(min(YM))),max(max(max(YM))),nE);
[XRM,YRM]=meshgrid(XR,YR);
% NaN-wrapping against the edge triangulation errors
XM=[2*XM(1,:,:)-XM(2,:,:); XM; 2*XM(end,:,:)-XM(end-1,:,:)]; XM=[2*XM(:,1,:)-XM(:,2,:) XM 2*XM(:,end,:)-XM(:,end-1,:)];
YM=[2*YM(1,:,:)-YM(2,:,:); YM; 2*YM(end,:,:)-YM(end-1,:,:)]; YM=[2*YM(:,1,:)-YM(:,2,:) YM 2*YM(:,end,:)-YM(:,end-1,:)];
Data=[NaN*Data(1,:,:); Data; NaN*Data(end,:,:)]; Data=[NaN*Data(:,1,:) Data NaN*Data(:,end,:)];

% re-mapping
DataR=zeros(nE,nA,nS);
for iS=1:nS
   iA=min(nA3,iS); X=XM(:,:,iA); X=X(:);  
   iE=min(nE3,iS); Y=YM(:,:,iE); Y=Y(:); 
   DL=Data(:,:,iS); DL=DL(:);
   % DataR(:,:,iS)=griddata(X,Y,D,AngleM,EnergyM);
   F=TriScatteredInterp(X,Y,DL); DataR(:,:,iS)=F(XRM,YRM);
end