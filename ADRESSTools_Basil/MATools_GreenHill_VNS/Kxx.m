function Kx = Kxx(HV,Eb,thtM,ThtA,surfNormX)
% Kx = Kxx(HV,Eb,thtM,ThtA,surfNormX) calculates k// along the analyzer slit 
% Inputs:
% HV (scalar or 1D vector) - photon energies
% Eb (Eb<0; scalar/column vector/2D/3D array) - binding energies
% thtM (scalar) - manipulator primary rotation
% ThtA (scalar/row vector/2D/3D array) - analyzer scale angles
% surfNormX - surface normal primary angle returned by SurfNormX.m
% IMPORTANT: The function implies the SX-ARPES@ADRESS standard geometry
% with the analyzer installed at 70deg to the incident X-rays

% parameters
global alpha ePhi thtMPhys nA nE
alpha=20; % nominal incidence angle
ePhi=4.5; % workfunction
thtMPhys=thtM+surfNormX;
nA=size(ThtA,2); nE=size(Eb,1); 
nS=max(length(HV),size(Eb,3));   % scan is either HV or tilt encoded in Eb

% conversion of iso-E 2D-arrays to 3D-arrays
if isscalar(Eb) && ~isvector(ThtA); ThtA=permute(ThtA,[3,2,1]); sq=1; else sq=0; end

% hv-independent data (scalar HV) including single image
Kx=[];
if isscalar(HV)
% - cycle through tilt
   for iS=1:nS
      EbF=Eb(:,:,min(iS,size(Eb,3))); ThtAF=ThtA;
      KxF=Kx1(HV,EbF,ThtAF);
      Kx=cat(3,Kx,KxF);
   end
% hv-dependent data (non-scalar HV)
else   
% - cycle through HV
   for iS=1:nS
      EbF=Eb(:,:,min(iS,size(Eb,3))); ThtAF=ThtA(:,:,min(iS,size(ThtA,3))); 
      KxF=Kx1(HV(iS),EbF,ThtAF);
      Kx=cat(3,Kx,KxF);
   end
end

% back conversion of iso-E 3D-arrays to 2D-arrays as required for ImData
if sq==1; Kx=(squeeze(Kx))'; end

% single image conversion
function Kx1 = Kx1(hv,Eb,ThtA)
% Kx1 = Kx1(hv,Eb,ThtA) calculates Kx for single Eb/ThtA image with scalar hv
global alpha ePhi thtMPhys nA nE
% Eb=squeeze(Eb); ThtA=squeeze(ThtA);
% forming 2D arrays unless scalars
if size(Eb,2)==1&&size(Eb,1)~=1; EbF=repmat(Eb,1,nA); else EbF=Eb; end
if size(ThtA,1)==1&&size(ThtA,2)~=1; ThtAF=repmat(ThtA,nE,1); else ThtAF=ThtA; end
% Kx calculation
Kx1=0.5124*sqrt(hv-ePhi+EbF).*sind(ThtAF+thtMPhys)-2*pi*hv*cosd(alpha+thtMPhys)/12400;