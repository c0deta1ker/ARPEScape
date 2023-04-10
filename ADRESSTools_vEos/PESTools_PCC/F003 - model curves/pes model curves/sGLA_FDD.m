function ydat = sGLA_FDD(xdat, x0, peak, fwhm, mr, asym, fdd_ef, fdd_T, fdd_fwhm)
% ydat = sGLA_FDD(xdat, x0, peak, fwhm, mr, asym, fdd_ef, fdd_T, fdd_fwhm)
%   Function that evaluates an asymmetric pseudo-Voigt curve with the
%   sGLA() function, after which it is multiplied by the Fermi-Dirac
%   Distribution FDDG(). These curve shapes are used for fitting EDCs near
%   the Fermi-edge, where the curve profiles may be chopped slightly due to
%   the close proximity of the Fermi-edge.
%
%   REQ. FUNCTIONS:
%   -   sGLA()
%   -   FDDG()
%
%   IN:
%   -   xdat:           N×1 column vector of the input domain (binding energy for PES).
%   -   x0:             scalar of the peak position along the x-axis of the Voigt.
%   -   peak:           scalar of the peak height of the Voigt.
%   -   fwhm:           scalar of the full-width at half-maximum (FWHM) of the Voigt.
%   -   mr:             mixing ratio; 0 for pure Gaussian, 1 is for pure Lorentzian
%   -   asym:           asymmetry ratio; 0 for none, 1 is for maximum
%   -   fdd_ef:         scalar of the Fermi-level position (eV).
%   -   fdd_T:          scalar of the temperature (K).
%   -   fdd_fwhm:       scalar of the full-width at half-maximum (FWHM) of the Gaussian broadening term (eV).
%
%   OUT:
%   -   ydat:           N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2;          x0 = 0; end
if nargin < 3;          peak = 1; end
if nargin < 4;          fwhm = 1; end
if nargin < 5;          mr = 0.5;  end
if nargin < 6;          asym = 0.5;  end
if nargin < 7;          fdd_ef = 0; end
if nargin < 8;          fdd_T = 12; end
if nargin < 9;          fdd_fwhm = 0; end
if isempty(x0);         x0 = 0; end
if isempty(peak);       peak = 1; end
if isempty(fwhm);       fwhm = 1; end
if isempty(mr);         mr = 0.5; end
if isempty(asym);       asym = 0.5; end
if isempty(fdd_ef);     fdd_ef = 0; end
if isempty(fdd_T);      fdd_T = 12; end
if isempty(fdd_fwhm);   fdd_fwhm = 0; end
%% Validity checks on the input parameters
if fwhm < 0; fwhm = 0; end
if mr < 0; mr = 0; end
if mr > 1; mr = 1; end
if asym < -1; asym = -1; end
if asym > 1; asym = 1; end
if fdd_T < 0; fdd_T = 0; end
if fdd_fwhm < 0; fdd_fwhm = 0; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% 1 - Determining the FDD and Gaussian functions to use
% - Extracting the intensity of the curve
int_sGLA  	= sGLA(xdat, x0, peak, fwhm, mr, asym);
% - Extracting the intensity of the FDD curve
int_FDDG 	= FDDG(xdat, fdd_ef, fdd_T, fdd_fwhm);

%% 2 - Multiplying the two functions together
ydat        = int_FDDG .* int_sGLA; 
ydat        = peak .* ydat / max(ydat);
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end