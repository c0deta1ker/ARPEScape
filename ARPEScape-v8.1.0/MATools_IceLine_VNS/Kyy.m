function Ky = Kyy(HV,Eb,ThtA,thtM,TltM,surfNormX,incAlpha,TltE)
% Ky = Kyy(HV,Eb,ThtA,thtM,TltM,surfNormX [,incAlpha] [,TltE]) calculates k// perpendicular to the analyzer slit
% Inputs:
% HV (scalar/vector) - photon energies
% Eb (Eb<0; scalar/column vector/2D/3D array) - binding energies
% ThtA (scalar/row vector/2D/3D array) - analyzer scale angles (may include angular warping)
% thtM - manipulator primary rotation
% TltM (scalar/row vector) - manipulator tilt relative to the surface normal
% surfNormX - surface normal primary angle returned by SurfNormX.m
% incAlpha - nominal X-ray incidence angle, for the present geometry equal to 9deg. If omitted, old geometry with 20deg implied.
% TltE (scalar/row vector) - angular deflection. If omitted, TltE=0. 
% Ver. 02.12.2022

% parameters
global alpha ePhi thtMPhys nA nE
if nargin<8; TltE=0; end % nominal incidence angle
if nargin<7; alpha=20; else alpha=incAlpha; end % nominal incidence angle
ePhi=4.5; % workfunction
nA=size(ThtA,2); nE=size(Eb,1); nH=length(HV); nT=length(TltM); nD=length(TltE);
if length(find([nH nT nD]>1))>1; disp('Error in Kyy: Only one of hv/tilt/defection can be vectorized'); return; end
if size(ThtA,1)==1&&size(ThtA,2)~=1; ThtA=repmat(ThtA,nE,1); end % form 2D array of ThtA unless scalar or 3D
if size(Eb,2)==1&&size(Eb,1)~=1; Eb=repmat(Eb,1,nA); end % form 2D array of Eb unless scalar or 3D
thtMPhys=thtM+surfNormX;
nS=max([length(HV) length(TltM) length(TltE)]); % scan is either HV or TltM or TltE

% conversion of iso-E 2D-arrays to 3D-arrays (immune to scalars)
% if isscalar(Eb); HV=permute(HV,[3,2,1]); TltM=permute(TltM,[3,2,1]); sq=1; else sq=0; end
if isscalar(Eb); ThtA=permute(ThtA,[3,2,1]); sq=1; else sq=0; end

% % hv-independent data (scalar HV) including single image
% Ky=[];
% if isscalar(HV)
% % - cycle through tilt (non-scalar TltM)
%    for iS=1:nS
%       KyF=Ky1(HV,Eb(:,:,min(iS,size(Eb,3))),ThtA,TltM(iS));
%       Ky=cat(3,Ky,KyF);
%    end
% % hv-dependent data (non-scalar HV)
% else   
% % - cycle through HV
%    for iS=1:nS
%       KyF=Ky1(HV(iS),Eb(:,:,min(iS,size(Eb,3))),ThtA(:,:,min(iS,size(ThtA,3))),TltM);
%       Ky=cat(3,Ky,KyF);
%    end
% end

% Ky conversion
Ky=[]; for iS=1:nS
      KyF=Ky1(HV(min(nH,iS)),squeeze(Eb(:,:,min(iS,size(Eb,3)))),squeeze(ThtA(:,:,min(iS,size(ThtA,3)))),TltM(min(nT,iS)),TltE(min(nD,iS)));      
% Note: by default all Matlab scalars and arrays have infinite number of singleton dimensions
   Ky=cat(3,Ky,KyF);
end
% back conversion of iso-E 3D-arrays to 2D-arrays as required for ImData
if sq==1; Ky=(squeeze(Ky))'; end

% single 2D-image conversion
function Ky1 = Ky1(hv,Eb,ThtA,tltM,tltE)
% Ky1 = Ky1(hv,Eb,ThtAF,tltM) calculates Ky for single Eb/ThtA image with scalar hv, tltM and tltE
global alpha ePhi thtMPhys nA nE
% forming 2D arrays if vectors
if size(Eb,2)==1&&size(Eb,1)~=1; Eb=repmat(Eb,1,nA); end
if size(ThtA,1)==1&&size(ThtA,2)~=1; ThtA=repmat(ThtA,nE,1); end
% Ky calculation
Ky1=0.5124*sqrt(hv-ePhi+Eb)*sind(tltM+tltE).*cosd(ThtA+thtMPhys)+2*pi*hv*sind(alpha+thtMPhys).*sind(tltM)/12400;