function int = Lorentzian(xdat, x0, peak, fwhm)
% int = Lorentzian(xdat, x0, peak, fwhm)
%   Function that evaluates a Lorentzian curve after defining its position,
%   full-width at half-maximum (FWHM) and amplitude.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       Nx1 column vector of the input domain (binding energy for XPS)
%   -   x0:         scalar of the peak position along the x-axis of the Lorentzian.
%   -   peak:       scalar of the total peak beneath the Lorentzian.
%   -   fwhm:     	scalar of the full-width at half-maximum (FWHM) of the Lorentzian.
%
%   OUT:
%   -   int:        Nx1 column vector of the output Lorentzian curve profile

%% Default parameters
if nargin < 4; fwhm = 1;  end
if nargin < 3; peak = 1; fwhm = 1;  end
if nargin < 2; x0 = 0; peak = 1; fwhm = 1; end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
if isempty(fwhm); fwhm = 1; end

%% - 1 - Determination of the Lorentzian curve intensities
% Ensuring xdat is a column vector
if size(xdat, 2) >1; xdat = xdat'; end
% Evaluating the normalised curve profile
int     = (peak .* 0.5*fwhm./pi) ./ ((xdat - x0).^2 + (0.5*fwhm).^2);
% Scaling the curve profile to match the desired peak
int 	= int ./ max(int(:));
int 	= peak .* int;
% If isnan, return zero
int(isnan(int)) = 0;
end