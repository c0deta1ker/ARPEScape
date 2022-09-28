function ECorr=CurveCorr(Angle,Energy,ep,CurveFit)
% ECorr=CurveCorr(Angle,Energy [,ep] [,CurveFit]) returns curvature corrected matrix
% of (binding or kinetic) energies ECorr for the vector of angles Angle, vector of 
% (binding or kinetic) energies Energy and pass energy ep with the curvature fit data 
% taken from the file CurveFit. If the input ep is omitted or empty, it is set to the 
% default value based on the input energy range.
% Note: CurveCorr assumes the raw size of the data array as produced by
% ReadARPES without trimming

disp('- Curvature correction')

% parameters
nA=length(Angle);
% default pass energy
if nargin<3; ep=[]; end
if isempty(ep)
   ep=(Energy(end)-Energy(1))/(2*0.066);
end
% default calibration file
if nargin<4; CurveFit=[]; end
if isempty(CurveFit)
   CurveFit='CurveFit.dat';
end

% mode WAM/MAM/LAD/MAD = 1/2/3/4
aRange=Angle(end)-Angle(1);
mode=4; if aRange>12; mode=3; end; if aRange>18; mode=2; end; if aRange>24; mode=1; end;

% load curvature correction
CurveFit=load(CurveFit);
% select the mode
ind=find(CurveFit(:,1)==mode); 
if isempty(ind); disp('- Error: Angular mode not curvature calibrated'); ECorr=[]; return
else CurveFit=CurveFit(ind,2:end);
end
% interpolate fitting coefficients
if size(CurveFit,1)==1
   curveFit=CurveFit(2:end); 
else
   Ep=[0;CurveFit(:,1);realmax];
   CurveFit=[CurveFit(1,2:end);CurveFit(:,2:end);CurveFit(end,2:end)];
   curveFit=interp1q(Ep',CurveFit,ep);
end
% append zero constant term
curveFit=[curveFit 0];

% angle-dependent energy correction
ECorr=polyval(curveFit,Angle)*ep;
% corrected energy matrix
ECorr=repmat(Energy,1,nA)-repmat(ECorr,length(Energy),1); 