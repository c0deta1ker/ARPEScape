function int = L2D(xdat, ydat, peak, x0, y0, fwhm)
% int = L2D(xdat, ydat, peak, x0, y0, fwhm)
%   Function that evaluates a 2D Lorentzian curve after defining its amplitude
%   and position in the x- and y-axis. The function assumes that the FWHM
%   is identical both in the x- and y-dimensions.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       [1×M] row vector of the x-axis input domain (kx for ARPES)
%   -   ydat:       [N×1] column vector of the y-axis input domain (Eb for ARPES)
%   -   peak:       scalar of the maximum peak height of the Lorentzian.
%   -   x0:         scalar of the peak position along the x-axis of the Lorentzian.
%   -   y0:         scalar of the peak position along the y-axis of the Lorentzian.
%   -   fwhm:       scalar of the full-width at half-maximum (FWHM) of the Lorentzian along both x- and y-axes.
%
%   OUT:
%   -   int:        [N×M] array of the output Lorentzian spot profile

%% Default parameters
if nargin < 6; fwhm = 1;  end
if nargin < 5; y0 = 0;  end
if nargin < 4; x0 = 0;  end
if nargin < 3; peak = 1;  end
if isempty(fwhm); fwhm = 1; end
if isempty(y0); y0 = 0; end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
% Validity check on inputs
if peak < 0; peak = 0; end
if fwhm < 0; fwhm = 0; end

%% - 1 - Determination of the Lorentzian curve intensities
% Ensuring xdat is a row vector
if size(xdat, 1) > 1 && size(xdat, 2) == 1; xdat = xdat'; end
% Ensuring ydat is a column vector
if size(ydat, 2) > 1 && size(ydat, 1) == 1; ydat = ydat'; end
% Evaluating the normalised curve profile
int     = (peak .* 0.5*fwhm./pi) ./ ((xdat - x0).^2 + (ydat - y0).^2 + (0.5*fwhm).^2).^1.5;
% Scaling the curve profile to match the desired peak
int 	= int ./ max(int(:));
int 	= peak .* int;
% If isnan, return zero
int(isnan(int)) = 0;
end