function ydat = PESCurve_FDD(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, FDEF, FDT, FDW, plot_result)
% ydat = PESCurve_FDD(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, FDEF, FDT, FDW, plot_result)
%   Function that evaluates a generic photoelectron spectroscopy (PES)
%   curve, by defining a primary (P) and spin-orbit split (SOS) component.
%   A Fermi-Dirac Distribution (FDD) is also multiplied with the PES curve.
%   This type of curve should be used to fit near Fermi-edge EDC cuts.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   cTYPE:      type of curve to use for fitting. Default: "sGLA" ("sGLA", "pGLA", "sGL", "pGL", "G", "L", "DS")
%   -   BE:      	scalar of the binding energy of PE curve.
%   -   INT:    	scalar of the peak intensity of PE curve.
%   -   FWHM:     	scalar of the FWHM of the PE curve.
%   -   MR:     	scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
%   -   LSE:     	scalar of the binding energy of spin-orbit split PE curve.
%   -   LSI:     	scalar of the branching ratio of spin-orbit split PE curve.
%   -   LSW:     	scalar of the additional lorentzian width of spin-orbit split PE curve.
%   -   ASY:     	scalar of the PE curve asymmetry parameter (usually for metallic systems).
%   -   FDEF:     	scalar of the FDD Fermi-Level position.
%   -   FDT:     	scalar of the FDD temperature.
%   -   FDW:     	scalar of the FDD Gaussian width after convolution.
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   ydat:     	N×1 column vector of the output FDD convolved Voigt curve profile

%% Default parameters
% Default based on inputs
if nargin < 2; cTYPE	= "sGLA"; end
if nargin < 3; BE   = 0.00; end
if nargin < 4; INT  = 1.00; end
if nargin < 5; FWHM = 0.25; end
if nargin < 6; MR = 0.25; end
if nargin < 7; LSE  = 0; end
if nargin < 8; LSI  = 0; end
if nargin < 9; LSW  = 0; end
if nargin < 10; ASY  = 0; end
if nargin < 11; FDEF  = 0; end
if nargin < 12; FDT  = 12; end
if nargin < 13; FDW  = 0.05; end
if nargin < 14; plot_result	= 0;  end
% Default based on empty inputs
if isempty(cTYPE);  cTYPE   = "sGLA"; end
if isempty(BE);     BE      = 0.00; end
if isempty(INT);    INT     = 1.00; end
if isempty(FWHM);   FWHM    = 0.25; end
if isempty(MR);     MR      = 0.25; end
if isempty(LSE);    LSE     = 0; end
if isempty(LSI);    LSI     = 0; end
if isempty(LSW);    LSW     = 0; end
if isempty(ASY);    ASY     = 0; end
if isempty(FDEF);   FDEF    = 0; end
if isempty(FDT);    FDT     = 12; end
if isempty(FDW);    FDW     = 0.05; end
if isempty(plot_result);    plot_result = 0; end
%% Validity checks on the input parameters
if INT < 0; INT = 0; end
if FWHM < 0; FWHM = 0; end
if MR < 0; MR = 0; end
if MR > 1; MR = 1; end
if LSI < 0; LSI = 0; end
if LSW < 0; LSW = 0; end
if ASY < -1; ASY = -1; end
if ASY > 1; ASY = 1; end
if FDW < 0; FDW = 0; end
if FDT < 0; FDT = 0; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% - 1 - Determination of the photoemission spectrum
% - Extracting the FDD curve
FDD_curve 	= FDDG(xdat, FDEF, FDT, FDW);
% - Extracting the PES curve
PES_int     = PESCurve(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY);
% - Multiplying the components together
ydat    	= PES_int .* FDD_curve;
% ydat         = INT .* (ydat / max(ydat));
% If isnan, return zero
ydat(isnan(ydat)) = 0;

%% - 2 - Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end
if size(ydat, 2) > 1; ydat = ydat'; end

%% -- For Debugging
if plot_result == 1
    % - Initialising the figure object
    figure(); hold on;
    % -- Plotting the 1D data
    plot(xdat, ydat, 'b-', 'linewidth', 2);
    plot(xdat, PES_int, 'k-', 'linewidth', 1);
    plot(xdat, FDD_curve, 'r-', 'linewidth', 2);
    gca_props(); title('PESCurve_FDD()', 'interpreter', 'none');  
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    legend({'PESCurve_FDD()', 'PESCurve()', 'FDDG()'}, 'location', 'best', 'Interpreter','none');
    % -- Determining the best limits for the plot
    Y = [ydat; PES_int; FDD_curve];
    axis([min(xdat(:)), max(xdat(:)), min(Y(:)), max(Y(:))]);
end

end