function fig = view_xps_init_model(xpsStr, cTYPE, iPESCurves, bTYPE, iPESBgrnd)
% fig = view_xps_init_model(xpsStr, cTYPE, iPESCurves, bTYPE, iPESBgrnd)
%   This function is used to plot the initial, model XPS curve PRIOR to
%   curve fitting with 'xps_solver()'. The plot consists of 3 subplots; (1) The
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
%   -   xpsStr:         MATLAB data-structure that contains the XPS data.
%   -   cTYPE:          1xN vector of the type of curve to use for the nth state.
%   -   iPESCurves:   	3 cells {x0}{lb}{ub} with Nx8 arrays: the n'th peak parameters [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%   -   iPESBgrnd:      3 cells {x0}{lb}{ub} with 1x5 vectors: the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%
%   OUT:
%   -   fig:  	MATLAB figure object with the ARPES data plotted.

%% - 1 - Extracting all background subtraction data and information
% -- Extracting the background
[X, D, B]  = PESBackground(xpsStr.xdat, xpsStr.int,...
    bTYPE, iPESBgrnd{1}(1), iPESBgrnd{1}(2), iPESBgrnd{1}(3), iPESBgrnd{1}(4), iPESBgrnd{1}(5), iPESBgrnd{1}(6));

%% - 2 - Extracting all model data and information
% -- Extracting the original XPS data
xdat   	= xpsStr.xdat;
int 	= xpsStr.int;
% -- Extracting the total number of states to be fitted
nSTATES = length(cTYPE);
% -- Initialising the vectors to contain XPS data
DB    = D - B;      	% Data after being background subtracted
M       = zeros(size(X));      	% XPS model data
% -- Filing through all curve components and extracting them seperately
cYY = [];
for i = 1:nSTATES
   	cYY(:,i) = PESCurve(X, cTYPE(i),...
        iPESCurves{1}(i,1), iPESCurves{1}(i,2), iPESCurves{1}(i,3),...
        iPESCurves{1}(i,4), iPESCurves{1}(i,5), iPESCurves{1}(i,6),...
        iPESCurves{1}(i,7), iPESCurves{1}(i,8));
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
    [iPESBgrnd{1}(1), iPESBgrnd{1}(1), iPESBgrnd{1}(2), iPESBgrnd{1}(2), iPESBgrnd{1}(1)],...
    [-1, 1, 1, -1, -1].*1e4, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', [0 0 0]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(xdat, int, 'b-', 'linewidth', 0.5*pp.lwidth);
plot(X, D, 'b-', 'linewidth', pp.llwidth);
plot(X, B, 'r-', 'linewidth', pp.llwidth);
plot(X, DB, 'k-', 'linewidth', pp.llwidth);
gca_props(); grid on;
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
axis([min(xdat(:)), max(xdat(:)),0, 1.25*max(int(:))]);
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
    line([iPESCurves{1}(i,1), iPESCurves{1}(i,1)], [0, max(cYY(:,i))],...
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

%% - 5 - Plotting the uncertainties to allow for easier optimisation
subplot(4,2,[2,4,6]); hold on;
for i = 1:nSTATES
    %% - 5.1 - Plotting the primary peak uncertainties
    lbBE  	= iPESCurves{2}(i,1); ubBE      = iPESCurves{3}(i,1);
    lbINT  	= iPESCurves{2}(i,2); ubINT     = iPESCurves{3}(i,2);
    x_vals = [lbBE, lbBE, ubBE, ubBE, lbBE];
    y_vals = [lbINT, ubINT, ubINT, lbINT, lbINT];
    patch(x_vals, y_vals, pp.col.fit{i}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);
    %% - 5.2 - Adding text to identify each peak index
    text(max(x_vals), max(y_vals), string(i), 'interpreter', 'latex', 'fontsize', 14, 'color', 'k');
    %% - 5.3 - Plotting the spin-orbit split (SOS) peak uncertainties
    BE      = iPESCurves{1}(i,1);   INT     = iPESCurves{1}(i,2);
    lbLSE  	= iPESCurves{2}(i,5);   ubLSE	= iPESCurves{3}(i,5);
    lbLSI  	= iPESCurves{2}(i,6);   ubLSI	= iPESCurves{3}(i,6);
    x_vals  = BE + [lbLSE, lbLSE, ubLSE, ubLSE, lbLSE];
    y_vals  = [lbINT*lbLSI, ubINT*ubLSI, ubINT*ubLSI, lbINT*lbLSI, lbINT*lbLSI];
    patch(x_vals, y_vals, pp.col.fit{i}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);
end
end