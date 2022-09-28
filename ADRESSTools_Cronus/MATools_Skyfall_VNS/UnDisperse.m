function [AA,ClipIm]=UnDisperse(A,hwSm,threshD2,threshD1)
% [AA,ClipIm]=UnDisperse(A,hwSm,threshD2,threshD1) clips from the ARPES image A the band-dispersion regions and
% evaluates the spectrum AA which is angle-integrated and normalized by the number of points outside these regions. 
% The dispersion regions are identified by the condition -d2I/dE2>threshD2 (between the peak inflexion points) or, 
% if the fourth parameter provided, by |d2I/dE2|>threshD2 and |dI/dE|>threshD1 (larger peak base width). Before 
% the derivation, the image A is Gaussian pre-smoothed in the energy direction with hwSm, and the resulting 
% spectrum AA post-smoothed with hwSm/2. If the -d2I/dE2>threshD2 condition is used, the corresponding peak 
% width underestimate can be compensted by an increase of hwSm (which simultaneously reduces noise).
% Example: [AA,ClipIm]=UnDisperse(Z,0.5/dE,0.04); ImData(X,Y,SetContrast(ClipIm,0.001,0.999)); plot(AA)

% parameters
histMax=0.999; Mask=ones(size(A));
% pre-smoothing
SmA=Gaco2(A,0,hwSm); 

% peak identification from second derivative only
if nargin<4
% - derivatives in the Y-direction
   DA2=-1*DifC(SmA,2); DA2(DA2<0)=0; DA2(isnan(DA2))=0;
% - normalization to histogrammic maximum
   DA2=DA2/max(max(DA2)); limDA2=stretchlim(DA2,[0 histMax]); DA2=DA2/limDA2(2);
% - clipping mask
   Mask(DA2>threshD2)=0;

% peak identification from second and first derivatives
else
% - derivatives in the Y-direction
   DA1=abs(DifC(SmA,1)); DA1(isnan(DA1))=0;
   DA2=-1*DifC(SmA,2); DA2(DA2<0)=0; DA2(isnan(DA2))=0; % peak width within -d2I/dY2 maximum 
%   DA2=abs(DifC(SmA,2)); DA2(isnan(DA2))=0; % peak width between -d2I/dY2 minima
% - normalization to histogrammic maximum
   DA1=DA1/max(max(DA1)); limDA1=stretchlim(DA1,[0 histMax]); DA1=DA1/limDA1(2);
   DA2=DA2/max(max(DA2)); limDA2=stretchlim(DA2,[0 histMax]); DA2=DA2/limDA2(2);
% - clipping mask
   Mask(DA2>threshD2|DA1>threshD1)=0;
end

% angle-integrated unstructure
% - input data with the clipped dispersion regions
ClipIm=A; ClipIm(Mask==0)=NaN;
% - nput data with zeros in the clipped regions  
A(Mask==0)=0;
% - angle-integrated data normalized by number of points outside the clipped regions   
nA=size(A,2); nE=size(A,1); AA=IntAngle(A,1:nA,(1:nE)',[1+0.05*(nA-1) 1+0.95*(nA-1)]);
% [~,AA]=Cut(1:nA,(1:nE)',A,'edc',[1+0.05*(nA-1) 1+0.95*(nA-1)]); % for some reasons does not work
MSum=IntAngle(Mask,1:nA,(1:nE)',[1+0.05*(nA-1) 1+0.95*(nA-1)]);
AA=AA./MSum;
% post-smoothing
AA=Gaco2(AA,0,hwSm/2);

% SmoZ=Gaco2(ExpZ,0,0.1/dE);
% figure; ImData(ExpX,ExpY,SmoZ); colormap hot; colorbar
% Di2Z=-DifC(SmoZ,2);
% Diff1=DifC(SmoZ,1);
% Cl2Z=Di2Z; Cl2Z(abs(Di2Z)>max(max(abs(Di2Z)))/40|abs(Diff1)>max(max(abs(Diff1)))/10)=NaN;
% figure; ImData(ExpX,ExpY,Cl2Z); colormap hot; colorbar;
% Cl2Z=Di2Z; Cl2Z(abs(Di2Z)>max(max(abs(Di2Z)))/20|abs(Diff1)>max(max(abs(Diff1)))/10)=NaN;
% figure; ImData(ExpX,ExpY,Cl2Z); colormap hot; colorbar;