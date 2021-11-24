function zdat = ARPESCurve(xdat, ydat, cTYPE, XLOC, YLOC, INT, FWHM, ASY)
% zdat = ARPESCurve(xdat, ydat, cTYPE, XLOC, YLOC, INT, FWHM, ASY)
%   Function that evaluates a generic angle-resolved photoelectron 
%   spectroscopy (ARPES) curve in 2D, as either a Gaussian or Lorentzian
%   function. These curves can be used to create simulated ARPES spectra.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	1xM row vector of the x-axis domain
%   -   ydat:    	Nx1 column vector of the y-axis domain
%   -   cTYPE:      type of curve to use for fitting. Default: "G2D" ("L2D", "G2DA")
%   -   XLOC:      	scalar of the x-location of the 2D PE curve.
%   -   YLOC:      	scalar of the y-location of the 2D PE curve.
%   -   INT:    	scalar of the peak intensity of 2D PE curve.
%   -   FWHM:     	scalar of the FWHM of the 2D PE curve.
%   -   ASY:     	scalar of the 2D PE curve asymmetry parameter (scales the FWHM along the y-axis).
%
%   OUT:
%   -   zdat:     	N x M matrix of the output 2D curve profile

%% Default parameters
% Default based on inputs
if nargin < 8; ASY  = 0;  end
if nargin < 7; FWHM = 0.25; end
if nargin < 6; INT  = 1.00; end
if nargin < 5; YLOC	= 0.00; end
if nargin < 4; XLOC	= 0.00; end
if nargin < 3; cTYPE = "G2D"; end
% Default based on empty inputs
if isempty(ASY);    ASY     = 0; end
if isempty(FWHM);   FWHM    = 0.25; end
if isempty(INT);    INT     = 1.00; end
if isempty(YLOC);   YLOC    = 0.00; end
if isempty(XLOC);   XLOC    = 0.00; end
if isempty(cTYPE);  cTYPE   = "G2D"; end
% Validity check on inputs
if INT < 0; INT = 0; end
if ASY < 0; ASY = 0; end
if FWHM < 0; FWHM = 0; end

%% - 1 - Determination of the photoemission spectrum
% Ensuring xdat is a row vector
if size(xdat, 1) > 1; xdat = xdat'; end
% Ensuring ydat is a column vector
if size(ydat, 2) > 1; ydat = ydat'; end
% Extracting the primary curve components
if cTYPE == "G2D"
    P_curve     = G2D(xdat, ydat, INT, XLOC, YLOC, FWHM); 
elseif cTYPE == "L2D"
    P_curve     = L2D(xdat, ydat, INT, XLOC, YLOC, FWHM); 
elseif cTYPE == "G2DA"
    P_curve     = G2DA(xdat, ydat, INT, XLOC, YLOC, FWHM, ASY.*FWHM); 
end
% Summing the components together
zdat = P_curve;
% If isnan, return zero
zdat(isnan(zdat)) = 0;

end