function int = pGL(xdat, x0, peak, fwhm, mr)
% int = pGL(xdat, x0, peak, fwhm, mr)
%   Function that evaluates a pseudo-Voigt curve by expressing the curve
%   profile as a PRODUCT of Gaussian/Lorentzian curve types.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	Nx1 column vector of the input domain (binding energy for XPS)
%   -   xc:         scalar of the peak position along the x-axis of the Voigt.
%   -   peak:    	scalar of the total peak beneath the Voigt.
%   -   fwhm:       scalar of the full-width at half-maximum (FWHM) of the Voigt.
%   -   mr:         mixing ratio; 0 for pure Gaussian, 1 is for pure Lorentzian
%
%   OUT:
%   -   int:    	Nx1 column vector of the output Voigt curve profile

%% Default parameters
if nargin < 5; mr = 0.5;  end
if nargin < 4; fwhm = 1; mr = 0.5; end
if nargin < 3; peak = 1; fwhm = 1; mr = 0.5;  end
if nargin < 2; x0 = 0; peak = 1; fwhm = 1; mr = 0.5;  end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
if isempty(fwhm); fwhm = 1; end
if isempty(mr); mr = 0.5; end

%% - 1 - Determination of the Voigt (pGL) curve intensities
% Ensuring xdat is a column vector
if size(xdat, 2) >1; xdat = xdat'; end
% Evaluating the normalised curve profile
int = exp(-4.*log(2).*(1-mr).*((xdat-x0).^2 ./ fwhm.^2)) ./ (1 + 4.*mr.*((xdat-x0).^2 ./ fwhm.^2));
% Scaling the curve profile to match the desired peak
int 	= int ./ max(int(:));
int 	= peak .* int;
% If isnan, return zero
int(isnan(int)) = 0;
end