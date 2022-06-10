function ydat = DoniachSunjic(xdat, x0, peak, fwhm, asym)
% ydat = DoniachSunjic(xdat, x0, peak, fwhm, asym)
%   Function that evaluates a Doniach-Sunjic curve after defining its position,
%   full-width at half-maximum (FWHM), amplitude and asymmetry. The formula
%   is used to characterise the asymmetry in the measured XPS of metallic
%   systems. The problem with this function is that it lacks proper
%   quantification, as F and E are both quantities that are related to the
%   FWHM and peak position respectively, but not exactly equal to them. It
%   is better to use Voigt lineshapes for fitting photoemission peaks.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   x0:         scalar of the (approximate) peak position.
%   -   peak:       scalar of the peak height.
%   -   fwhm:       scalar of the (approximate) FWHM.
%   -   asym:       asymmetry ratio; 0 for none, 1 is for maximum
%
%   OUT:
%   -   ydat:       N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2; x0 = 0; end
if nargin < 3; peak = 1; end
if nargin < 4; fwhm = 1; end
if nargin < 5; asym = 0.5;  end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
if isempty(fwhm); fwhm = 1; end
if isempty(asym); asym = 0.5; end
%% Validity checks on the input parameters
% -- If the FWHM is negative, set it to zero
if fwhm < 0; fwhm = 0; end
% -- If the asym is negative, set it to zero
if asym < 0; asym = 0; end
% -- If the asym is positive, set it to one
if asym > 1; asym = 1; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% - 1 - Determination of the Doniach-Sunjic curve intensities
% Evaluating the normalised curve profile
ydat     = (cos(0.5*pi*asym + (1-asym).*atan((xdat - x0) / fwhm))) ./ ((fwhm.^2+(xdat-x0).^2).^((1-asym)/2));
% Scaling the curve profile to match the desired peak
ydat 	= ydat ./ max(ydat(:));
ydat 	= peak .* ydat;
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end