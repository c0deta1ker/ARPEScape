function ydat = sGL(xdat, x0, peak, fwhm, mr)
% ydat = sGL(xdat, x0, peak, fwhm, mr)
%   Function that evaluates a pseudo-Voigt curve by expressing the curve
%   profile as a SUM of Gaussian/Lorentzian curve types.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   x0:         scalar of the peak position along the x-axis of the Voigt.
%   -   peak:       scalar of the peak height of the Voigt.
%   -   fwhm:       scalar of the full-width at half-maximum (FWHM) of the Voigt.
%   -   mr:         mixing ratio; 0 for pure Gaussian, 1 is for pure Lorentzian
%
%   OUT:
%   -   ydat:       N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2;      x0 = 0; end
if nargin < 3;      peak = 1; end
if nargin < 4;      fwhm = 1; end
if nargin < 5;      mr = 0.5;  end
if isempty(x0);     x0 = 0; end
if isempty(peak);   peak = 1; end
if isempty(fwhm);   fwhm = 1; end
if isempty(mr);     mr = 0.5; end
%% Validity checks on the input parameters
% -- If the FWHM is negative, set it to zero
if fwhm < 0; fwhm = 0; end
% -- If the MR is negative, set it to zero
if mr < 0; mr = 0; end
% -- If the MR is > 1, set it to 1
if mr > 1; mr = 1; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% - 1 - Determination of the Voigt (pGL) curve intensities
% Evaluating the normalised curve profile
ydat    = (1-mr).*exp(-4.*log(2).*((xdat-x0).^2 ./ fwhm.^2)) + mr ./ (1 + 4.*((xdat-x0).^2 ./ fwhm.^2));
% Scaling the curve profile to match the desired peak
ydat 	= ydat ./ max(ydat(:));
ydat 	= peak .* ydat;
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end