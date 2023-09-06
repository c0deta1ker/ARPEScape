function ydat = PESCurve(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, plot_result)
% ydat = PESCurve(xdat, cTYPE, BE, INT, FWHM, MR, LSE, LSI, LSW, ASY, plot_result)
%   Function that evaluates a generic photoelectron spectroscopy (PES)
%   curve, by defining a primary (P) and spin-orbit split (SOS) component.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	N×1 column vector of the input domain (binding energy for XPS)
%   -   cTYPE:      type of curve to use for fitting. Default: "sGLA" ("sGLA", "pGLA", "sGL", "pGL", "G", "L", "DS")
%   -   BE:      	scalar of the binding energy of PE curve.
%   -   INT:    	scalar of the peak intensity of PE curve.
%   -   FWHM:     	scalar of the FWHM of the PE curve.
%   -   MR:     	scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
%   -   LSE:     	scalar of the binding energy of spin-orbit split PE curve.
%   -   LSI:     	scalar of the branching ratio of spin-orbit split PE curve.
%   -   LSW:     	scalar of the additional lorentzian width of spin-orbit split PE curve.
%   -   ASY:     	scalar of the PE curve asymmetry parameter (usually for metallic systems).
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   ydat:     	N×1 column vector of the output Voigt curve profile

%% Default parameters
% Default based on inputs
if nargin < 2;  cTYPE	= "sGLA"; end
if nargin < 3;  BE      = 0.00; end
if nargin < 4;  INT     = 1.00; end
if nargin < 5;  FWHM    = 0.25; end
if nargin < 6;  MR      = 0.50; end
if nargin < 7;  LSE     = 0;  end
if nargin < 8;  LSI     = 0;  end
if nargin < 9;  LSW     = 0;  end
if nargin < 10; ASY     = 0;  end
if nargin < 11; plot_result	= 0;  end
% Default based on empty inputs
if isempty(cTYPE);  cTYPE   = "sGLA"; end
if isempty(BE);     BE      = 0.00; end
if isempty(INT);    INT     = 1.00; end
if isempty(FWHM);   FWHM    = 0.25; end
if isempty(MR);     MR      = 0.50; end
if isempty(LSE);    LSE     = 0; end
if isempty(LSI);    LSI     = 0; end
if isempty(LSW);    LSW     = 0; end
if isempty(ASY);    ASY     = 0; end
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

%% - 1 - Determination of the PES curve
% - Extracting the primary and spin-split components
if strcmpi(cTYPE,"sGLA")
    P_curve     = sGLA(xdat, BE, INT, FWHM, MR, ASY);
    SOS_curve 	= sGLA(xdat, BE+LSE, LSI*INT, FWHM+LSW, MR, ASY);    
elseif strcmpi(cTYPE,"pGLA")
    P_curve     = pGLA(xdat, BE, INT, FWHM, MR, ASY);
    SOS_curve 	= pGLA(xdat, BE+LSE, LSI*INT, FWHM+LSW, MR, ASY);  
elseif strcmpi(cTYPE,"sGL")
    P_curve     = sGL(xdat, BE, INT, FWHM, MR);
    SOS_curve 	= sGL(xdat, BE+LSE, LSI*INT, FWHM+LSW, MR);    
elseif strcmpi(cTYPE,"pGL")
    P_curve     = pGL(xdat, BE, INT, FWHM, MR);
    SOS_curve 	= pGL(xdat, BE+LSE, LSI*INT, FWHM+LSW, MR);
elseif strcmpi(cTYPE,"G")
    P_curve     = Gauss(xdat, BE, INT, FWHM);
    SOS_curve 	= Gauss(xdat, BE+LSE, LSI*INT, FWHM+LSW);    
elseif strcmpi(cTYPE,"L")
    P_curve     = Lorentzian(xdat, BE, INT, FWHM);
    SOS_curve 	= Lorentzian(xdat, BE+LSE, LSI*INT, FWHM+LSW);
elseif strcmpi(cTYPE,"DS")
    P_curve     = DoniachSunjic(xdat, BE, INT, FWHM, ASY);
    SOS_curve 	= DoniachSunjic(xdat, BE+LSE, LSI*INT, FWHM+LSW, ASY);
else
    error('cTYPE: Curve type not recognised.');
end
% - Summing the components together
ydat = P_curve + SOS_curve;
% - If isnan, return zero
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
    gca_props(); title('PESCurve()', 'interpreter', 'none'); 
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    legend({'PESCurve()'}, 'location', 'best', 'fontsize', 9);
    % -- Determining the best limits for the plot
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
end
end