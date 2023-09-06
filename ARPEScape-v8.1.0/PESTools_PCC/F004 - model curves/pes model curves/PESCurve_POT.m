function [ydat, ydat2D] = PESCurve_POT(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, MFP, ZPOT, EPOT, plot_result)
% [ydat, ydat2D] = PESCurve_POT(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, MFP, ZPOT, EPOT, plot_result)
%   Function that evaluates a generic photoelectron spectroscopy (PES)
%   curve, by defining a primary (P) and spin-orbit split (SOS) component.
%   A potential profile is also mapped onto the curve
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
%   -   MFP:     	scalar of the mean-free path of the emitted photoelectrons (either from optical, TPP-2M or fits).
%   -   ZPOT:     	1×M array of the z-domain (depth) of the potential profile.
%   -   EPOT:     	1×M array of the potential energy relative to the Fermi-level.
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   ydat:     	N×1 column vector of the output PES curve profile.
%   -   ydat2D:    	N×M array of all the PES curve profiles at each point in the potential profile.

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
if nargin < 11; MFP  = 1; end
if nargin < 12; EPOT  = 0; end
if nargin < 13; ZPOT  = 0; end
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
if isempty(MFP);    MFP     = 1; end
if isempty(EPOT);   EPOT    = 0; end
if isempty(ZPOT);   ZPOT  	= 0; end
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
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% - 1 - Determination of the photoemission spectrum
% - Setting the potential to be zero at the BE value
EPOT    = BE + EPOT;  % EPOT    = BE + (EPOT - EPOT(1));
% - Extracting the PES curve
for i = 1:length(ZPOT)
    ydat2D(:,i) = PESCurve(xdat, cTYPE, EPOT(i), INT.*exp(-ZPOT(i)./MFP), FWHM, MR, LSE, LSI, LSW, ASY);
end
% - Normalising all the individual curves for comparison
ydat2D     = INT .* (ydat2D ./ max(ydat2D(:)));
% - Extracting the final intensity
ydat         = sum(ydat2D, 2);
ydat         = INT .* (ydat / max(ydat));
% If isnan, return zero
ydat(isnan(ydat)) = 0;

%% - 2 - Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end
if size(ydat, 2) > 1; ydat = ydat'; end

%% -- For Debugging
% -- Initialising the figure
if plot_result == 1
    fig = figure(); set(fig, 'position', [10, 10, 900, 400]);
    cols = jet(length(ZPOT)+2);
    % -- Plotting all individual curves, shifted by the potential and the final sum
    subplot(121); hold on;
    % --- Plotting individual curves
    for i = 1:length(ZPOT)
        plot(xdat, ydat2D(:,i), 'k-', 'color', cols(i,:), 'linewidth', 0.5);
    end
    % --- Plotting the final curve
    plot(xdat, ydat, 'k-', 'linewidth', 2);
    % --- Formatting the figure
    gca_props(); title('PESCurve_POT()', 'Interpreter','none');
    xlabel(' E_B - E_F (eV) ', 'fontweight', 'bold');
    ylabel(' Intensity ', 'fontweight', 'bold');
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
    % -- Plotting the potential profile information
    subplot(122); hold on;
    % --- Image of the curve series shifted by potential and scaled by MFP
    ImData(ZPOT, xdat, ydat2D);
    % --- Plotting the potential energy curves
    plot(ZPOT, EPOT, 'r-', 'linewidth', 2);
    plot(ZPOT, BE*ones(size(ZPOT)), 'g-', 'linewidth', 1);
    % --- Formatting the figure
    img_props(); cbar_props('jet');
    xlabel(' z (nm) ', 'fontweight', 'bold');
    ylabel(' E_B - E_F (eV) ', 'fontweight', 'bold');
    axis([0, max(ZPOT(:)), min(xdat(:)), max(xdat(:))]);
end
end