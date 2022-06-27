function ydat = pGLA(xdat, x0, peak, fwhm, mr, asym)
% ydat = pGLA(xdat, x0, peak, fwhm, mr, asym)
%   Function that evaluates a, asymmetric pseudo-Voigt curve by expressing 
%   the curve profile as a PRODUCT of Gaussian/Lorentzian curve types.
%   The curve is modified here with the asymmetric, exponential blend, 
%   where:  pGLA(x) = pGL(x) + AEB(x) * (1 - pGL(x)). 
%   AEB(x) is determined from (AEB()).
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   x0:         scalar of the peak position along the x-axis of the Voigt.
%   -   peak:       scalar of the peak height of the Voigt.
%   -   fwhm:       scalar of the full-width at half-maximum (FWHM) of the Voigt.
%   -   mr:         mixing ratio; 0 for pure Gaussian, 1 is for pure Lorentzian
%   -   asym:       asymmetry ratio; 0 for none, 1 is for maximum
%
%   OUT:
%   -   ydat:       N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2;      x0 = 0; end
if nargin < 3;      peak = 1; end
if nargin < 4;      fwhm = 1; end
if nargin < 5;      mr = 0.5;  end
if nargin < 6;      asym = 0.5;  end
if isempty(x0);     x0 = 0; end
if isempty(peak);   peak = 1; end
if isempty(fwhm);   fwhm = 1; end
if isempty(mr);     mr = 0.5; end
if isempty(asym);   asym = 0.5; end
%% Validity checks on the input parameters
if fwhm < 0; fwhm = 0; end
if mr < 0; mr = 0; end
if mr > 1; mr = 1; end
if asym < -1; asym = -1; end
if asym > 1; asym = 1; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% - 1 - Determination of the Voigt (pGL) curve intensities
int_pGL     = pGL(xdat, x0, peak, fwhm, mr);
%% - 2 - Determination of the Exponential Asymmetric Blend (AEB) intensities
int_AEB     = AEB(xdat, x0, asym);
%% - 3 - Determination of pGLA curve intensities
ydat        = int_pGL + int_AEB .* (max(int_pGL(:)) - int_pGL);
% Scaling the curve profile to match the desired peak
ydat        = ydat ./ max(ydat(:));
ydat        = peak .* ydat;
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end