function ydat = Lorentzian(xdat, x0, peak, fwhm)
% ydat = Lorentzian(xdat, x0, peak, fwhm)
%   Function that evaluates a Lorentzian curve after defining its position,
%   full-width at half-maximum (FWHM) and amplitude.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   x0:         scalar of the peak position along the x-axis of the Lorentzian.
%   -   peak:       scalar of the peak height of the Lorentzian.
%   -   fwhm:       scalar of the full-width at half-maximum (FWHM) of the Lorentzian.
%
%   OUT:
%   -   ydat:       N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2; x0 = 0; end
if nargin < 3; peak = 1; end
if nargin < 4; fwhm = 1; end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
if isempty(fwhm); fwhm = 1; end
%% Validity checks on the input parameters
% -- If the FWHM is negative, set it to zero
if fwhm < 0; fwhm = 0; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% - 1 - Determination of the Lorentzian curve intensities
% Evaluating the normalised curve profile
ydat    = (peak .* 0.5*fwhm./pi) ./ ((xdat - x0).^2 + (0.5*fwhm).^2);
% Scaling the curve profile to match the desired peak
ydat 	= ydat ./ max(ydat(:));
ydat 	= peak .* ydat;
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end