function ydat = Gauss(xdat, x0, peak, fwhm)
% ydat = Gauss(xdat, x0, peak, fwhm)
%   Function that evaluates a Gaussian curve after defining its position,
%   full-width at half-maximum (FWHM) and amplitude.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   x0:         scalar of the peak position along the x-axis of the Gaussian.
%   -   peak:       scalar of the peak height of the Gaussian.
%   -   fwhm:       scalar of the full-width at half-maximum (FWHM) of the Gaussian.
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

%% - 1 - Determination of the Gaussian curve intensities
% Using the standard deviation as the input width so that the FWHM is correct
sigma   = fwhm ./ (2*sqrt(2*log(2))); 
% Evaluating the normalised curve profile
ydat    = (1 ./ (sigma * sqrt(2*pi))) .* exp(-0.5.*((xdat-x0)./sigma).^2);
% Scaling the curve profile to match the desired peak
ydat 	= ydat ./ max(ydat(:));
ydat 	= peak .* ydat;
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end