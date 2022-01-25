function [int, PES_int] = PESCurve_POT(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, MFP, ZPOT, EPOT)
%int = PESCurve_POT(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, MFP, ZPOT, EPOT)
%   Function that evaluates a generic photoelectron spectroscopy (PES)
%   curve, by defining a primary (P) and spin-orbit split (SOS) component.
%   A potential profile is also mapped onto the curve
%   This type of curve should be used to fit near Fermi-edge EDC cuts.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	Nx1 column vector of the input domain (binding energy for XPS)
%   -   cTYPE:      type of curve to use for fitting. Default: "sGLA" ("pGLA", "sGL", "pGL", "G", "L", "DS")
%   -   BE:      	scalar of the binding energy of PE curve.
%   -   INT:    	scalar of the peak intensity of PE curve.
%   -   FWHM:     	scalar of the FWHM of the PE curve.
%   -   MR:     	scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
%   -   LSE:     	scalar of the binding energy of spin-orbit split PE curve.
%   -   LSI:     	scalar of the branching ratio of spin-orbit split PE curve.
%   -   LSW:     	scalar of the additional lorentzian width of spin-orbit split PE curve.
%   -   ASY:     	scalar of the PE curve asymmetry parameter (usually for metallic systems).
%   -   MFP:     	scalar of the mean-free path of the emitted photoelectrons (either from optical, TPP-2M or fits).
%   -   ZPOT:     	1xN array of the z-domain (depth) of the potential profile.
%   -   EPOT:     	1xN array of the potential energy relative to the Fermi-level.
%
%   OUT:
%   -   int:     	Nx1 column vector of the output FDD convolved Voigt curve profile

%% Default parameters
% Default based on inputs
if nargin < 13; ZPOT  = 0;  end
if nargin < 12; EPOT  = 0;  end
if nargin < 11; MFP  = 1;  end
if nargin < 10; ASY  = 0;  end
if nargin < 9; LSW  = 0;  end
if nargin < 8; LSI  = 0;  end
if nargin < 7; LSE  = 0;  end
if nargin < 6; MR = 0.25; end
if nargin < 5; FWHM = 0.25; end
if nargin < 4; INT  = 1.00; end
if nargin < 3; BE   = 0.00; end
if nargin < 2; cTYPE	= "sGLA"; end
% Default based on empty inputs
if isempty(ZPOT);   ZPOT  	= 0; end
if isempty(EPOT);   EPOT    = 0; end
if isempty(MFP);    MFP     = 1; end
if isempty(ASY);    ASY     = 0; end
if isempty(LSW);    LSW     = 0; end
if isempty(LSI);    LSI     = 0; end
if isempty(LSE);    LSE     = 0; end
if isempty(MR);     MR      = 0.25; end
if isempty(FWHM);   FWHM    = 0.25; end
if isempty(INT);    INT     = 1.00; end
if isempty(BE);     BE      = 0.00; end
if isempty(cTYPE);  cTYPE   = "sGLA"; end

%% - 1 - Determination of the photoemission spectrum
% Ensuring xdat is a column vector
if size(xdat, 2) >1; xdat = xdat'; end
% Validity check on inputs
if INT < 0; INT = 0; end
if ASY < 0; ASY = 0; end
if FWHM < 0; FWHM = 0; end
% - Setting the potential to be zero at the BE value
EPOT    = BE + (EPOT - EPOT(1));
% - Extracting the PES curve
for i = 1:length(ZPOT)
    PES_int(:,i) = PESCurve(xdat, cTYPE, EPOT(i), INT.*exp(-ZPOT(i)./MFP), FWHM, MR, LSE, LSI, LSW, ASY);
end
% - Extracting the final intensity
int         = sum(PES_int, 2);
int         = INT .* (int / max(int));
% If isnan, return zero
int(isnan(int)) = 0;

% %% - 2 - Plotting the results for debugging
% % -- Initialising the figure
% fig = figure(); set(fig, 'position', [10, 10, 900, 400]);
% cols = jet(length(ZPOT)+2);
% % -- Plotting all individual curves, shifted by the potential and the final sum
% subplot(121); hold on;
% % --- Plotting individual curves
% for i = 1:length(ZPOT)
%     plot(xdat, PES_int(:,i), 'k-', 'color', cols(i,:), 'linewidth', 0.5);
% end
% % --- Plotting the final curve
% plot(xdat, int, 'k-', 'linewidth', 2);
% % --- Formatting the figure
% gca_props(); title('PESCurve_POT()', 'Interpreter','none');
% xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
% ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
% axis([min(xdat(:)), max(xdat(:)), min(int(:)), max(int(:))]);
% % -- Plotting the potential profile information
% subplot(122); hold on;
% % --- Image of the curve series shifted by potential and scaled by MFP
% ImData(ZPOT, xdat, PES_int);
% % --- Plotting the potential energy curves
% plot(ZPOT, EPOT, 'r-', 'linewidth', 2);
% plot(ZPOT, BE*ones(size(ZPOT)), 'g-', 'linewidth', 1);
% % --- Formatting the figure
% img_props(); colormap jet;
% xlabel('$$ \bf  z (nm) $$', 'Interpreter', 'latex');
% ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
% axis([0, max(ZPOT(:)), min(xdat(:)), max(xdat(:))]);

end