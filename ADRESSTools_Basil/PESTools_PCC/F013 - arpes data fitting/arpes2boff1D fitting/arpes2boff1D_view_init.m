function fig = arpes2boff1D_view_init(arpesStr, cTYPE, FUNC, iparams, bTYPE, ibgrnd)
% fig = arpes2boff1D_view_init(arpesStr, cTYPE, FUNC, iparams, bTYPE, ibgrnd)
%   This function is used to plot the initial, model curves PRIOR to
%   curve fitting with 'arpes2boff1D_solver()'. The plot consists of 3 subplots; (1) The
%   background that is determined from the fit; (2) A plot showing all of
%   the fitted curve components, as well as the final model fit and
%   experimental data; (3) A plot of the residuals, showing the quality of
%   the experimental and model fit. This is used as an informative plot
%   that allows you to view and create a better initial guess of the XPS
%   model prior to running the fitting algorithm.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   arpesStr:       MATLAB data-structure that contains the initial ARPES data.
%   -   FUNC:           1xN cell of functions that gives the band offset vs subband energy
%   -   cTYPE:          1xN vector of the type of curve to use for fitting. Default: "sGLA" ("pGLA", "DS")
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x6 array: the model fit parameters [FDEF,FDT,FDW,BOFF,INT,FWHM]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x5 vectors: the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%
%   OUT:
%   -   fig:  	MATLAB figure object with the ARPES data plotted.

%% Initialising variables
% -- Extracting the total number of states to be fitted
nSTATES = length(cTYPE);
% -- Extracting the input variables for the fit
FDEF    = iparams{1}(1);        % scalar of the FDD Fermi-Level position.
FDT     = iparams{1}(2);        % scalar of the FDD temperature.
FDW     = iparams{1}(3);        % scalar of the FDD Gaussian width after convolution.
BOFF    = iparams{1}(4);        % scalar of the band offset
MR      = iparams{1}(5);        % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
INT     = iparams{1}(6:6+nSTATES-1);             % scalar of the peak intensity of PE curve.
FWHM    = iparams{1}(6+nSTATES:6+2*nSTATES-1);	 % scalar of the FWHM of the PE curve.

%% - 1 - Extracting all background subtraction data and information
% -- Extracting the background
[X, D, B]  = PESBackground(arpesStr.xdat, arpesStr.ydat,...
    bTYPE, ibgrnd{1}(1), ibgrnd{1}(2), ibgrnd{1}(3), ibgrnd{1}(4), ibgrnd{1}(5), ibgrnd{1}(6));

%% - 2 - Extracting all model data and information
% -- Extracting the original XPS data
xdat   	= arpesStr.xdat;
ydat 	= arpesStr.ydat;
% -- Initialising the vectors to contain XPS data
DB    = D - B;                  % Data after being background subtracted
M     = zeros(size(X));      	% XPS model data
% -- Filing through all curve components and extracting them seperately
cYY = []; BE = [];
for i = 1:nSTATES
    BE(i) = FUNC{i}(BOFF);
    cYY(:,i) = PESCurve_FDD(X, cTYPE(i), BE(i), INT(i), FWHM(i), MR, 0, 0, 0, 0, FDEF, FDT, FDW);
    M = M + cYY(:,i);
end

%% - 3 - Determination of the residuals and chi-squared
R           = M - (D - B);	% Residuals
CHISQ       = sum(R.^2 ./ abs(M));          % Chi-squared

%% - 4 - Plotting the model to be used for fitting
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'XPS Curve Fitting');
fig.Position(3) = 2.5*pp.fig5x4(1); 
fig.Position(4) = pp.fig5x4(2);
%% - 4.1 - PLOTTING THE RAW DATA AND BEST FIT BACKGROUND
subplot(121); hold on;
h = patch(...
    [ibgrnd{1}(1), ibgrnd{1}(1), ibgrnd{1}(2), ibgrnd{1}(2), ibgrnd{1}(1)],...
    [-1, 1, 1, -1, -1].*1e4, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', [0 0 0]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(xdat, ydat, 'b-', 'linewidth', 0.5*pp.lwidth);
plot(X, D, 'b-', 'linewidth', pp.llwidth);
plot(X, B, 'r-', 'linewidth', pp.llwidth);
plot(X, DB, 'k-', 'linewidth', pp.llwidth);
gca_props(); grid on;
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
axis([min(xdat(:)), max(xdat(:)),0, 1.25*max(ydat(:))]);
legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
title("Background Subtraction", 'interpreter', 'none', 'fontsize', 9);
%% - 4.2 - PLOTTING THE BEST FIT CURVE COMPONENTS AND FINAL MODEL FIT
subplot(4,2,[2,4,6]); hold on;
% -- Plotting all of the curve components
for i = 1:nSTATES
    area(X, cYY(:,i), 'FaceColor', pp.col.fit{i}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
    plot(X, cYY(:,i), 'k-', 'linewidth', 0.25);
end
% -- Plotting the curve component energy locations
for i = 1:nSTATES
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    line([1, 1]*BE(i), [0, max(cYY(:,i))],...
        'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '-');
end
% -- Plotting the experimental and fit spectra
plot(X, DB, 'k-', 'color', pp.col.dat{1}, 'linewidth', 2*pp.llwidth);
plot(X, M, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Add annotation for the quality of fit
text(0.04, 0.90, "$$ \chi^2 = $$ " + string(CHISQ),...
    'interpreter', 'latex', 'fontsize', 14, 'color', 'k', 'Units','normalized');
% - Formatting the figure
gca_props(); grid on;
ylabel('$$ \bf  Intensity$$', 'Interpreter', 'latex');
axis([min(X(:)), max(X(:)), min(DB(:)), 1.10*max(DB(:))]);
title("Initial Curves for Fitting", 'interpreter', 'none', 'fontsize', 9);
%% - 4.3 - PLOTTING THE RESIDUALS FOR THE BEST FIT TO THE DATA
subplot(4,2,8); hold on;
bar(X, R, 'facecolor', [0 0 0]);
gca_props(); grid on;
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');

end