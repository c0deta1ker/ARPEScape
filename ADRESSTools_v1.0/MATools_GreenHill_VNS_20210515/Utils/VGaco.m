function AA=VGaco(A,HW,hs,PadMethod)
% AA=VGaco(A,HW [,hs] [,PadMethod]) convolutes A (vector) with Gaussian 
% having varying half-width HW (vector or scalar, in inter-point distances) 
% and fixed half-support hs (scalar, in inter-point distances). If skipped 
% or empty, the default hs=3*hw. Reduction of the end effects is controlled 
% by PadMethod (optional): default symmetric padding if skipped or 'symm', 
% or padding by hs-wide arrays linearly extrapolated from the HW-wide end 
% intervals if 'extrap'

% default support
if nargin<3; hs=[]; end
if isempty(hs); hs=3*HW; end

% check the PadMethod arguments
if nargin==4 && ~strcmp(PadMethod,'extrap') && ~strcmp(PadMethod,'symm') 
   error('''symm'' (default) and ''extrap'' padding methods are supported'); 
end

% parameters
na=length(A); AA=A; A=reshape(A,[1 na]);
HW=HW/sqrt(log(2))+eps; if length(HW)==1; HW=HW*ones(1,na); end
hs=ceil(hs);

% convolution with padding by linear extrapolation
if nargin==4 && strcmp(PadMethod,'extrap')
% - linear fit within HW
   rangeL=1:1+ceil(HW(1)); PL=polyfit(rangeL,A(rangeL),1); AL=polyval(PL,1-hs:0);
   rangeR=na-ceil(HW(end)):na; PR=polyfit(rangeR,A(rangeR:na),1); AR=polyval(PR,na+1:na+hs);
% - appending
   if size(A,1)<size(A,2) A=[AL A AR]; else A=[AL; A; AR]; end
% - convolution
   for i=1:na
      RANGE=i-hs:i+hs;
      WEIGHTS=exp( -1*((i-RANGE)/HW(i)).^2 );
      AA(i)=dot( A(RANGE+hs),WEIGHTS )/sum( WEIGHTS ); 
   end

% convolution with symmetric padding (default)
else
% - convolution with re-weighting at the ends
   for i=1:na
      RANGE=max([1 i-hs]):min([na i+hs]);
      WEIGHTS=exp( -1*((i-RANGE)/HW(i)).^2 );
      AA(i)=dot( A(RANGE),WEIGHTS )/sum( WEIGHTS ); 
   end 
end

% end of VGaco