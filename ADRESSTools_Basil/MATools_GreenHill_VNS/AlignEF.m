function [EAlign,EF,Fail]=AlignEF(Data,ECorr,eWin,dEWin,dESmooth,feat)
% [EAlign,EF,Fail]=AlignEF(Data,ECorr [,eWin] [,dEWin] [,dESmooth] [,feat]) aligns the energy scale relative 
% to EF or peak . The inputs are 2D or 3D data array Data, corresponding 1D-column % vector or 
% (curvature corrected) 2D-matrix of energies ECorr, (optional scalar/vector of the length size(Data,3) or 
% empty value) approximate EF position eWin, (optional or empty value) width of the EF support dEWin, 
% (optional or empty value) Gaussian pre-smoothing energy width dESmooth, and (optional or empty 
% value) feat = 'Peak' (case-insensitive) defines the feature to align as a peak instead of the default Fermi 
% edge. If eWin,dEWin and dESmooth are skipped or empty values, they take reasonable automatic 
% values.  The output EF is the determined EF positions in the input energy scale. Fail returns the indices 
% of the scans where the alignment failed and EF was set to eWin.
% - Ver. 7 Mar 2020 (added vectorization of eWin)
% disp('- EF alignment')

% parameters
% - dimensions
nE=size(Data,1); nA=size(Data,2); nS=size(Data,3);
% - energy matrix
if size(ECorr,2)==1; rep=1; ECorr=repmat(ECorr,1,nA); else rep=0; end
% - energy interval and step
eMin=max(ECorr(1,:)); eMax=min(ECorr(end,:)); eInterval=eMax-eMin; dE=eInterval/(size(ECorr,1)-1);
% - failed alignment indices
Fail=[];
% automatic parameters
if nargin<6; feat='Edge'; end
if nargin<5; dESmooth=[]; end; if isempty(dESmooth); dESmooth=0.01*eInterval; end
if nargin<4; dEWin=[]; end; if isempty(dEWin); dEWin=0.1*eInterval; end
if nargin<3; eWin=[]; end; if isempty(eWin); eWin=eMin+0.7*eInterval; dEWin=0.1*eInterval; end
if isscalar(eWin); eWin=eWin*ones(1,size(Data,3)); end
if length(eWin)~=size(Data,3); disp('- AlignEF: Incompatible eWin length'); return; end

% angle integration (can be replaced by evalc(IntAngle(...), where evalc suppresses text output)
% - reference grid
EInt=(eMin:dE:eMax)';
% - cycle over angles in triads for speed
AngleInt=zeros(length(EInt),nS);
for iS=1:nS
   if nA>1    
      for iA=round(0.05*nA)+1:3:round(0.95*nA)-1
% - interpolation on the reference grid       
         ACorr=interp1(ECorr(:,iA),sum(Data(:,iA-1:iA+1,iS),2),EInt);
% - angle integration      
         AngleInt(:,iS)=AngleInt(:,iS)+ACorr;
      end
   else
      AngleInt(:,iS)=Data(:,1,nS);
   end
end

% processing
% - logarithmic derivative if edge not peak
if ~isequal(lower(feat),'peak')
  D=AngleInt(3:length(EInt),:)-AngleInt(1:length(EInt)-2,:); D=[D(1,:);D;D(end,:)];
  D=D./AngleInt;
else D=-1*AngleInt; end
% - smoothing (after the derivation to suppress the spikes)
hwPts=dESmooth/dE; hsPts=3*hwPts; DS=Gaco1(D,hwPts,hsPts);

% alignment at EF
% - cycle through scans
EAlign=zeros(nE,nA,nS); EF=zeros(1,nS); 
for iS=1:nS
   % - alignment region 
   DS1=DS(EInt>=eWin(iS)-dEWin/2&EInt<=eWin(iS)+dEWin/2,iS);
   ERef=EInt(EInt>=eWin(iS)-dEWin/2&EInt<=eWin(iS)+dEWin/2);
% - locating the minimum
   L=DS1(1:end-2); C=DS1(2:end-1); R=DS1(3:end); IMinima=1+find(C<L&C<R);
   if isempty(IMinima)
      eF=eWin(iS); Fail=[Fail iS];
      disp(['Error in AlignEF.m: Failed to locate EF, scan = ' num2str(iS)]);
   else
% - eF as the global minimum 
      [~,iMin]=min(DS1(IMinima)); ind=IMinima(iMin);
      x=XYMax(-1*(DS1(ind-1:ind+1))');
      eF=ERef(ind-1)+x*(ERef(ind)-ERef(ind-1));
   end
% - energy matrix alignment
   EAlign(:,:,iS)=ECorr-eF; EF(iS)=eF;
end

% - reduction to energy vector
if rep==1; EAlign=mean(EAlign,2); end