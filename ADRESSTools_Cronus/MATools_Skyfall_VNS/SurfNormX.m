function surfNormX = SurfNormX(hv,eB,kx,thtM,thtA,incAlpha)
% surfNormX = SurfNormX(hv,eB,kx,thtM,thtA,incAlpha) calculates the surface normal
% angle in the measurement plane. Inputs: 
% hv - photon energy;
% eB and kx - binding energy (eB<0) and known k// of some reference feature;
% thtN and thtA - corresponding angle of the sample rotation and angle on 
% the analyser scale.
% incAlpha - nominal X-ray incidence angle, for the present geometry equal to 9deg. If omitted, old geometry with 20deg implied.
% Ver. 08.11.2022

% parameters
if nargin<6; alpha=20; else alpha=incAlpha; end % nominal incidence angle
ePhi=4.5; % workfunction
% inputs check
if isempty(hv*eB*kx*thtM*thtA); surfNormX=[]; return; end
% physical reference THT
[thtMPhys,~,exitFlag]=fsolve(@(x) 0.5124*sqrt(hv-ePhi+eB)*sind(thtA+x)-...
                    2*pi*hv*cosd(alpha+x)/12400-kx,0,optimset('Display','off','TolFun',1e-10));
if exitFlag~=1||imag(thtMPhys)~=0; thtMPhys=[]; msgbox('Error: Failed to find the surface normal'); end
% surface normal
surfNormX=thtMPhys-thtM;