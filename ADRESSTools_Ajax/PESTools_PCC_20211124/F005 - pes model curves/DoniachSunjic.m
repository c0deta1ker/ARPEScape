function int = DoniachSunjic(xdat, x0, peak, fwhm, asym)
% int = DoniachSunjic(xdat, x0, peak, fwhm, asym)
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
%   -   xdat:    	Nx1 column vector of the input domain (binding energy for XPS)
%   -   x0:         scalar of the (approximate) peak position.
%   -   peak:    	scalar of the total peak beneath the Voigt.
%   -   fwhm:       scalar of the (approximate) FWHM.
%   -   asym:       asymmetry ratio; 0 for none, 1 is for maximum
%
%   OUT:
%   -   int:        Nx1 column vector of the output Gaussian curve profile

%% Default parameters
if nargin < 5; asym = 0.5;  end
if nargin < 4; fwhm = 1; asym = 0.5; end
if nargin < 3; peak = 1; fwhm = 1; asym = 0.5;  end
if nargin < 2; x0 = 0; peak = 1; fwhm = 1; asym = 0.5;  end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
if isempty(fwhm); fwhm = 1; end
if isempty(asym); asym = 0.5; end

%% - 1 - Determination of the Doniach-Sunjic curve intensities
% Ensuring xdat is a column vector
if size(xdat, 2) > 1; xdat = xdat'; end
% Evaluating the normalised curve profile
int     = (cos(0.5*pi*asym + (1-asym).*atan((xdat - x0) / fwhm))) ./ ((fwhm.^2+(xdat-x0).^2).^((1-asym)/2));
% Scaling the curve profile to match the desired peak
int 	= int ./ max(int(:));
int 	= peak .* int;
% If isnan, return zero
int(isnan(int)) = 0;
end