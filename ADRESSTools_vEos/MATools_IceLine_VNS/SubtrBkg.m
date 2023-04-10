function A=SubtrBkg(ACorr,EAlign,Data,EWinNorm,scCoeff,clipNeg)
% A=SubtrBkg(ACorr,EAlign,Data [,EWinNorm] [,scCoeff] [,clipNeg]) subtracts angle-integrated spectral weight from
% the 2D/3D data array Data. Before the subtraction, the data is normalized to the spectral intensity integrated within
% the energy window EWinNorm (default 90% of the EAlign range). scCoeff scales the subtracted angle-integrated
% intensity (default 0.7) and clipNeg indicates whether the negative values are set to zero (default 1) or not (0).  
% Ver. 11.07.2022

% default arguments
 if nargin<6||isempty(clipNeg); clipNeg=1; end
 if nargin<5||isempty(scCoeff); scCoeff=0.7; end
 if nargin<4||isempty(EWinNorm); EWinNorm=[]; end
 
 % validity check
if ~isnumeric(scCoeff)||~isequal(size(scCoeff),[1 1]); disp ('Error in SubtrBkg: scCoeff should be scalar'); A=[]; return; end
if ~isequal(clipNeg,0)&&~isequal(clipNeg,1); disp ('Error in SubtrBkg: clipNeg should be 0 or 1'); A=[]; return; end
if ~isempty(EWinNorm)&&~isequal(size(EWinNorm),[1 2])&&~isequal(size(EWinNorm),[2 1]); disp('Error in SubtrBkg: EWinNorm should be a two-element vector'); A=[]; return; end

% normalization energy window
ERange=[min(min(EAlign(:,1,:))) max(max(EAlign(:,end,:)))];
if isempty(EWinNorm)
   EWinNorm=[0.9*ERange(1)+0.1*ERange(2) 0.1*ERange(1)+0.9*ERange(2)];
else
   if min(EWinNorm)<ERange(1)||max(EWinNorm)>ERange(2) disp('Error in SubtrBkg: Inconsistent normalization energy window'); A=[]; return; end
end
% normalizing data
disp('- Normalizing data')
IDiv=Slice(ACorr,EAlign,Data,'isoE',EWinNorm);
IDiv=Gaco2(IDiv,5,0); % HWHM adjusted to suppress uneven sensitivity of the CCD channels
IDiv=repmat(IDiv,1,1,size(Data,1)); IDiv=permute(IDiv,[3 2 1]);
Data=Data./IDiv;

 % angle-integration window
aMin=max(max(ACorr(:,1,:))); aMax=min(min(ACorr(:,end,:))); 
AWin=[0.95*aMin+0.05*aMax 0.05*aMin+0.95*aMax];
% subtract angle-integrated spectrum and clip negative values
if scCoeff~=0
   disp('- Subtracting angle-integrated background')
   ISubtr=Slice(ACorr,EAlign,Data,'isoK',AWin);
   ISubtr=repmat(ISubtr,1,1,size(Data,2)); ISubtr=permute(ISubtr,[1 3 2]);
   A=Data-scCoeff*ISubtr;
   if clipNeg==1; A(A<0)=0; end
end