function int = sGLA_FDD(xdat, x0, peak, fwhm, mr, asym, fdd_ef, fdd_T, fdd_fwhm)
% int = sGLA_FDD(xdat, x0, peak, fwhm, mr, asym, fdd_ef, fdd_T, fdd_fwhm)
%   Function that evaluates an asymmetric pseudo-Voigt curve with the
%   sGLA() function, after which it is multiplied by the Fermi-Dirac
%   Distribution FDDG(). These curve shapes are used for fitting EDCs near
%   the Fermi-edge, where the curve profiles may be chopped slightly due to
%   the close proximity of the Fermi-edge.
%
%   REQ. FUNCTIONS:
%   -   sGLA(xdat, x0, peak, fwhm, mr, asym)
%   -   FDDG(xdat, ef, T, fwhm)
%
%   IN:
%   -   xdat:    	Nx1 column vector of the input domain (binding energy for XPS)
%   -   x0:         scalar of the peak position along the x-axis of the Voigt.
%   -   peak:    	scalar of the total peak beneath the Voigt.
%   -   fwhm:       scalar of the full-width at half-maximum (FWHM) of the Voigt.
%   -   mr:         mixing ratio; 0 for pure Gaussian, 1 is for pure Lorentzian
%   -   asym:       asymmetry ratio; 0 for none, 1 is for maximum
%   -   fdd_ef:   	scalar of the Fermi-level position (eV).
%   -   fdd_T:    	scalar of the temperature (K).
%   -   fdd_fwhm: 	scalar of the full-width at half-maximum (FWHM) of the Gaussian broadening term (eV).
%
%   OUT:
%   -   int:        Nx1 column vector of the output Voigt curve profile

%% Default parameters
if nargin < 9; fdd_fwhm = 0;  end
if nargin < 8; fdd_T = 12;  end
if nargin < 7; fdd_ef = 0;  end
if nargin < 6; asym = 0.1;  end
if nargin < 5; mr = 0.5; asym = 0.1; end
if nargin < 4; fwhm = 1; mr = 0.5; asym = 0.1; end
if nargin < 3; peak = 1; fwhm = 1; mr = 0.5; asym = 0.1; end
if nargin < 2; x0 = 0; peak = 1; fwhm = 1; mr = 0.5; asym = 0.1; end
if isempty(x0); x0 = 0; end
if isempty(peak); peak = 1; end
if isempty(fwhm); fwhm = 1; end
if isempty(mr); mr = 0.5; end
if isempty(asym); asym = 0.1; end
if isempty(fdd_ef); fdd_ef = 0; end
if isempty(fdd_T); fdd_T = 12; end
if isempty(fdd_fwhm); fdd_fwhm = 0; end

%% 1 - Determining the FDD and Gaussian functions to use
% Ensuring xdat is a column vector
if size(xdat, 2) >1; xdat = xdat'; end
% - Extracting the intensity of the curve
int_sGLA  	= sGLA(xdat, x0, peak, fwhm, mr, asym);
% - Extracting the intensity of the FDD curve
int_FDDG 	= FDDG(xdat, fdd_ef, fdd_T, fdd_fwhm);

%% 2 - Multiplying the two functions together
int      = int_FDDG .* int_sGLA; 
int      = peak .* int / max(int);
% If isnan, return zero
int(isnan(int)) = 0;
end