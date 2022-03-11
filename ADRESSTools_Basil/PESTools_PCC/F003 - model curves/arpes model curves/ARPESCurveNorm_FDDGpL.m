function zdat = ARPESCurveNorm_FDDGpL(xdat, ydat, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR, FDEF, FDT, FDW, BGR, BIN, BCO)
% zdat = ARPESCurve_FDD(xdat, ydat, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR, FDEF, FDT, FDW)
%   Function that evaluates a generic angle-resolved photoelectron 
%   spectroscopy (ARPES) curve in 2D. This can be used to model an ARPES
%   spectrum, that is built up by a sum of individual Gaussians.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	1xM row vector of the x-axis domain
%   -   ydat:    	Nx1 column vector of the y-axis domain
%   -   cTYPE:      type of curve to use for fitting. Default: "G2DA" ("G2D", "L2D")
%   -   INT:    	scalar of the peak intensity of 2D PE curve.
%   -   XLOC:      	scalar of the x-location of the min/max of the 2D parabolic ARPES dispersion [Ang^-1].
%   -   YLOC:      	scalar of the y-location of the min/max of the 2D parabolic ARPES dispersion [eV].
%   -   XFWHM:     	scalar of the x-axis FWHM for each Gaussian (k-resolution) [Ang^-1]
%   -   YFWHM:     	scalar of the y-axis FWHM for each Gaussian (Eb-resolution) [eV]
%   -   MSTAR:     	scalar of the effective mass, which governs the curvature of the parabola.
%   -   FDEF:     	scalar of the FDD Fermi-Level position.
%   -   FDT:     	scalar of the FDD temperature.
%   -   FDW:     	scalar of the FDD Gaussian width after convolution.
%   -   BGR:        scalar of the gradient of the linear background.
%   -   BIN:        scalar of the y-intercept of the linear background.
%   -   BCO:        scalar of the constant background y-offset value.
%
%   OUT:
%   -   zdat:     	N x M matrix of the output 2D parabolic ARPES dispersion.

%% Default parameters
% Default based on inputs
if nargin < 15; BCO  = 0.00;  end
if nargin < 14; BIN  = 1.00;  end
if nargin < 13; BGR  = 0.00;  end
if nargin < 12; FDW  = 0.05;  end
if nargin < 11; FDT  = 12;  end
if nargin < 10; FDEF  = 0;  end
if nargin < 9; MSTAR = 0.20;  end
if nargin < 8; YFWHM = 0.25; end
if nargin < 7; XFWHM = 0.02; end
if nargin < 6; YLOC	= 0.00; end
if nargin < 5; XLOC	= 0.00; end
if nargin < 4; INT  = 1.00; end
if nargin < 3; cTYPE = "G2DA"; end
% Default based on empty inputs
if isempty(cTYPE);  cTYPE   = "G2DA"; end
if isempty(INT);    INT     = 1.00; end
if isempty(XLOC);   XLOC    = 0.00; end
if isempty(YLOC);   YLOC    = 0.00; end
if isempty(XFWHM);  XFWHM   = 0.02; end
if isempty(YFWHM);  YFWHM   = 0.25; end
if isempty(MSTAR);  MSTAR   = 0.20; end
if isempty(FDW);    FDW     = 0.05; end
if isempty(FDT);    FDT     = 12; end
if isempty(FDEF);   FDEF    = 0; end
if isempty(BGR);    BGR     = 0; end
if isempty(BIN);    BIN     = 1; end
if isempty(BCO);    BCO     = 0; end

%% - 1 - Determination of the angle-resolved photoemission spectrum
% Ensuring xdat is a row vector
if size(xdat, 1) > 1; xdat = xdat'; end
% Ensuring ydat is a column vector
if size(ydat, 2) > 1; ydat = ydat'; end
% - Extracting the FDD curve
FDD_curve 	= FDDG(ydat, FDEF, FDT, FDW);
% - Extracting the ARPES curve
ARPES_curve	= ARPESCurveNorm(xdat, ydat, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR);
% - Multiplying the components together
zdat        = (ARPES_curve + BGR.*ydat+BIN) .* FDD_curve + BCO;
% int         = INT .* (int / max(int));
% If isnan, return zero
zdat(isnan(zdat)) = 0;

end