function zdat = ARPESCurveNorm_FDD(xdat, ydat, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR, FDEF, FDT, plot_result)
% zdat = ARPESCurveNorm_FDD(xdat, ydat, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR, FDEF, FDT, plot_result)
%   Function that evaluates a generic angle-resolved photoelectron 
%   spectroscopy (ARPES) curve in 2D. This can be used to model an ARPES
%   spectrum, that is built up by a sum of individual Gaussians. The
%   intensity is normalised across all EDCs.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       [1×M] row vector of the x-axis input domain (kx for ARPES)
%   -   ydat:       [N×1] column vector of the y-axis input domain (Eb for ARPES)
%   -   cTYPE:      type of curve to use for fitting. Default: "G2DA" ("G2D", "L2D")
%   -   INT:    	scalar of the peak intensity of 2D PE curve.
%   -   XLOC:      	scalar of the x-location of the min/max of the 2D parabolic ARPES dispersion [Ang^-1].
%   -   YLOC:      	scalar of the y-location of the min/max of the 2D parabolic ARPES dispersion [eV].
%   -   XFWHM:     	scalar of the x-axis FWHM for each Gaussian (k-resolution) [Ang^-1]
%   -   YFWHM:     	scalar of the y-axis FWHM for each Gaussian (Eb-resolution) [eV]
%   -   MSTAR:     	scalar of the effective mass, which governs the curvature of the parabola.
%   -   FDEF:     	scalar of the FDD Fermi-Level position.
%   -   FDT:     	scalar of the FDD temperature.
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   zdat:     	[N×M] matrix of the output 2D parabolic ARPES dispersion.

%% Default parameters
% Default based on inputs
if nargin < 12; plot_result = 0;  end
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
if isempty(FDT);    FDT     = 12; end
if isempty(FDEF);   FDEF    = 0; end
if isempty(plot_result);  plot_result   = 0; end

%% - 1 - Determination of the angle-resolved photoemission spectrum
% Ensuring xdat is a row vector
if size(xdat, 1) > 1; xdat = xdat'; end
% Ensuring ydat is a column vector
if size(ydat, 2) > 1; ydat = ydat'; end
% - Extracting the FDD curve
FDD_curve 	= FDD(ydat, FDEF, FDT);
% - Extracting the ARPES curve
ARPES_curve	= ARPESCurveNorm(xdat, ydat, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR);
% - Multiplying the components together
zdat        = ARPES_curve .* FDD_curve;
% int         = INT .* (int / max(int));
% If isnan, return zero
zdat(isnan(zdat)) = 0;

%% -- For Debugging
if plot_result == 1
    pp = plot_props();
    fig = figure(); 
    fig.Position(3) = pp.fig4x4(1); fig.Position(4) = pp.fig4x4(2);
    hold on;
    ImData(xdat, ydat, zdat); 
    img_props(); cbar_props(); title('ARPESCurveNorm_FDD()', 'Interpreter', 'none');
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
end

end