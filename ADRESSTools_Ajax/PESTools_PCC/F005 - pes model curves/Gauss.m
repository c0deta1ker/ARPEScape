function int = Gauss(xdat, x0, peak, fwhm)
% int = Gauss(xdat, x0, peak, fwhm)
%   Function that evaluates a Gaussian curve after defining its position,
%   full-width at half-maximum (FWHM) and amplitude.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       Nx1 column vector of the input domain (binding energy for XPS)
%   -   x0:         scalar of the peak position along the x-axis of the Gaussian.
%   -   peak:       scalar of the peak height the Gaussian.
%   -   fwhm:     	scalar of the full-width at half-maximum (FWHM) of the Gaussian.
%
%   OUT:
%   -   int:        Nx1 column vector of the output Gaussian curve profile

%% Default parameters
if nargin < 4; fwhm = 1;  end
if nargin < 3; peak = 1; fwhm = 1;  end
if nargin < 2; x0 = 0; peak = 1; fwhm = 1; end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
if isempty(fwhm); fwhm = 1; end

%% - 1 - Determination of the Gaussian curve intensities
% Ensuring xdat is a column vector
if size(xdat, 2) > 1; xdat = xdat'; end
% Using the standard deviation as the input width so that the FWHM is correct
sigma   = fwhm ./ (2*sqrt(2*log(2))); 
% Evaluating the normalised curve profile
int     = (1 ./ (sigma * sqrt(2*pi))) .* exp(-0.5.*((xdat-x0)./sigma).^2);
% Scaling the curve profile to match the desired peak
int 	= int ./ max(int(:));
int 	= peak .* int;
% If isnan, return zero
int(isnan(int)) = 0;
end