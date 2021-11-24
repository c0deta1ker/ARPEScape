function int = pGLA(xdat, x0, peak, fwhm, mr, asym)
% int = pGLA(xdat, x0, peak, fwhm, mr, asym)
%   Function that evaluates a, asymmetric pseudo-Voigt curve by expressing 
%   the curve profile as a SUM of Gaussian/Lorentzian curve types.
%   The curve is modified here with the asymmetric, exponential blend, 
%   where:  pGLA(x) = pGL(x) + AEB(x) * (1 - pGL(x)). 
%   AEB(x) is determined from (AEB()).
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	Nx1 column vector of the input domain (binding energy for XPS)
%   -   x0:         scalar of the peak position along the x-axis of the Voigt.
%   -   peak:    	scalar of the total peak beneath the Voigt.
%   -   fwhm:       scalar of the full-width at half-maximum (FWHM) of the Voigt.
%   -   mr:         mixing ratio; 0 for pure Gaussian, 1 is for pure Lorentzian
%   -   asym:       asymmetry ratio; 0 for none, 1 is for maximum
%
%   OUT:
%   -   int:        Nx1 column vector of the output Voigt curve profile

%% Default parameters
if nargin < 6; asym = 0.1;  end
if nargin < 5; mr = 0.5; asym = 0.1; end
if nargin < 4; fwhm = 1; mr = 0.5; asym = 0.1; end
if nargin < 3; peak = 1; fwhm = 1; mr = 0.5; asym = 0.1; end
if nargin < 2; x0 = 0; peak = 1; fwhm = 1; mr = 0.5; asym = 0.1; end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
if isempty(fwhm); fwhm = 1; end
if isempty(mr); mr = 0.5; end
if isempty(asym); asym = 0.1; end

%% - 1 - Determination of the Voigt (pGL) curve intensities
int_sGL     = pGL(xdat, x0, peak, fwhm, mr);
%% - 2 - Determination of the Exponential Asymmetric Blend (T) intensities
int_AEB     = AEB(xdat, x0, asym);
%% - 3 - Determination of pGLA curve intensities
int         = int_sGL + int_AEB .* (max(int_sGL(:)) - int_sGL);
% Scaling the curve profile to match the desired peak
int 	= int ./ max(int(:));
int 	= peak .* int;
% If isnan, return zero
int(isnan(int)) = 0;
end