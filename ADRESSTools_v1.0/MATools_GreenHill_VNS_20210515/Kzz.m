function Kz = Kzz(HV,Eb,thtM,ThtA,TltM,v000,surfNormX)
% Kz = Kzz(HV,Eb,thtM,ThtA,TltM,v000,surfNormX) calculates kz along the surface normal 
% Inputs:
% HV (scalar or vector) - photon energies
% Eb (Eb<0; scalar/column vector/2D/3D array) - binding energies
% ThtA (scalar/row vector/2D/3D array) - analyzer scale angles
% TltM (scalar/row vector) - manipulator tilt relative to its normal emission value  
% ThtA and TltM (scalar/column vector/2D/3D array) - angles along the x- and y-axis
% v000 >0 - inner potential relative to EF
% thtM - the manipulator primary rotation
% surfNormX - surface normal primary angle returned by SurfNormX.m
% IMPORTANT: The function implies the SX-ARPES@ADRESS standard geometry and 
% axes notation with the analyser slit oriented in the MP.
% Ver. 02-05-2019

% parameters
global alpha ePhi thtMPhys v0001 nA nE
alpha=20; % nominal incidence angle
ePhi=4.5; % workfunction
thtMPhys=thtM+surfNormX;
v0001=abs(v000);
nA=size(ThtA,2); nE=size(Eb,1); nH=length(HV); nT=length(TltM);
nS=max(nH,size(TltM,3));   % number of scans in HV or in tilt encoded in TltM
if size(ThtA,1)==1&&size(ThtA,2)~=1; ThtA=repmat(ThtA,nE,1); end % form 2D array of ThtA unless scalar or 3D
if size(Eb,2)==1&&size(Eb,1)~=1; Eb=repmat(Eb,1,nA); end % form 2D array of Eb unless scalar or 3D

% conversion of iso-E/iso-K 2D-arrays from Slice to 3D-arrays (immune to scalars)
if isscalar(Eb); ThtA=permute(ThtA,[3,2,1]); TltM=permute(TltM,[3,2,1]); sq=1; else sq=0; end
if isscalar(ThtA); Eb=permute(Eb,[1,3,2]); TltM=permute(TltM,[1,3,2]); sq=2; end

% % extension of Eb to 3D (if not extended by AlignEF)
% if nS>1 && size(Eb,3)==1; Eb=repmat(Eb,1,1,nS); end
% % extension of ThtA,TltM to 3D
% if nS>1 && size(ThtA,3)==1; ThtA=repmat(ThtA,1,1,nS); end
% if nS>1 && size(TltM,3)==1; TltM=repmat(TltM,1,1,nS); end

% Kz conversion
Kz=[]; for iS=1:nS
   KzF=Kz1(HV(min(nH,iS)),squeeze(Eb(:,:,min(iS,size(Eb,3)))),squeeze(ThtA(:,:,min(iS,size(ThtA,3)))),TltM(min(nT,iS)));  
% Note: by default all Matlab scalars and arrays have infinite number of singleton dimensions
   Kz=cat(3,Kz,KzF);
end

% back conversion of iso-E/iso-K 3D-arrays to 2D-arrays as required for ImData
if sq==1; Kz=(squeeze(Kz))'; end
if sq==2; Kz=squeeze(Kz); end

% single 2D-image conversion
function Kz1 = Kz1(hv,Eb,ThtA,tltM)
% Kz1 = Kz1(hv,Eb,ThtA,tltM) calculates Kz for single Eb/ThtA image with scalar hv
global alpha ePhi thtMPhys v0001 nA nE
% forming 2D arrays if vectors
if size(Eb,2)==1&&size(Eb,1)~=1; Eb=repmat(Eb,1,nA); end
if size(ThtA,1)==1&&size(ThtA,2)~=1; ThtA=repmat(ThtA,nE,1); end
% final-state Kz
Kx1_2=(hv-ePhi+Eb).*sind(ThtA+thtMPhys).^2;
Ky1_2=(hv-ePhi+Eb).*(sind(tltM).*cosd(ThtA+thtMPhys)).^2; 
Kz1=0.5124*sqrt(hv+Eb+v0001-(Kx1_2+Ky1_2));
% photon momentum corrected initial-state Kz
Kz1=Kz1+2*pi*hv*sind(alpha+thtMPhys)/12400;