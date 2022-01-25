function fig = kz2model_view_fit(fitStr)
% fig = kz2model_view_fit(fitStr)
%   This function is used to plot the initial, kz line profile model PRIOR to
%   curve fitting with 'kz2model_solver()'. The plot consists of 2 subplots;
%   (1) A plot showing all of the fitted curve components, as well as the 
%   final model fit and experimental data; (2) A plot of the residuals, 
%   showing the quality of the experimental and model fit. This is used as an 
%   informative plot that allows you to view and create a better initial 
%   guess of the model prior to running the fitting algorithm.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   xdat:           1xN vector that contains the domain of the line profile.
%   -   ydat:           1xN vector that contains the intensity data of the line profile.
%   -   cTYPE:          string of the type of curves to use; "sGLA" ("pGLA", "sGL", "pGL", "G", "L", "DS").
%   -   iparams:        3 cells {x0}{lb}{ub} with Nx8 arrays: the n'th peak parameters [KZ0,INT,MFP,DZ,MR,SKZ0,SINT,ASY,BGR]
%
%   OUT:
%   -   fig:  	MATLAB figure object with the ARPES data plotted.

%% - 1 - Extracting all model data and information
X   = fitStr.xdat;     	% Original domain
D   = fitStr.ydat;     	% Original data
M0  = fitStr.M0;      	% Model data (MFP)
M   = fitStr.M;         % Model data (MFP+KZ)
B   = fitStr.B;         % Background
DB	= D - B;           	% Data after being background subtracted
XX  = fitStr.XX;        % Interpolated domain
YY0	= fitStr.YY0;       % Interpolated MFP broadened spectrum
YY1	= fitStr.YY1;       % Interpolated Kz broadened spectrum
YY  = fitStr.YY;        % Interpolated MFP + kz broadened spectrum
%% - 3 - Determination of the residuals and chi-squared
R           = M - (D - B);             % Residuals
CHISQ       = sum(R.^2 ./ abs(M));   	% Chi-squared

%% - 4 - Plotting the model to be used for fitting
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'XPS Curve Fitting');
fig.Position(3) = 1.5*pp.fig5x4(1); 
fig.Position(4) = pp.fig5x4(2);
%% - 4.2 - PLOTTING THE BEST FIT CURVE COMPONENTS AND FINAL MODEL FIT
subplot(4,1,[1,2,3]); hold on;
% -- Plotting all of the curve components
area(XX, YY, 'FaceColor', pp.col.fit{2}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
plot(XX, YY, 'k-', 'linewidth', 0.25);
area(XX, YY1, 'FaceColor', pp.col.fit{4}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
plot(XX, YY1, 'k-', 'linewidth', 0.25);
area(XX, YY0, 'FaceColor', pp.col.fit{3}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
plot(XX, YY0, 'k-', 'linewidth', 0.25);
% -- Plotting the curve component energy locations
line([fitStr.KZ0, fitStr.KZ0], [0, fitStr.INT], 'Color', [0 0 0], 'LineWidth', 2, 'Linestyle', '-');
line(fitStr.KZ0+[fitStr.SKZ0, fitStr.SKZ0], [0, fitStr.INT*fitStr.SINT], 'Color', [0 0 0], 'LineWidth', 2, 'Linestyle', '-');
% -- Plotting the experimental and fit spectra
plot(X, DB, 'k-', 'color', pp.col.dat{1}, 'linewidth', 2*pp.llwidth);
plot(XX, YY, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Add annotation for the quality of fit
text(0.04, 0.90, "$$ \chi^2 = $$ " + string(CHISQ),...
    'interpreter', 'latex', 'fontsize', 14, 'color', 'k', 'Units','normalized');
% - Formatting the figure
gca_props(); grid on;
ylabel('$$ \bf  Intensity$$', 'Interpreter', 'latex');
axis([min(X(:)), max(X(:)), min(DB(:)), 1.10*max(DB(:))]);
title("Best Fit Outcome", 'interpreter', 'none', 'fontsize', 9);
%% - 4.3 - PLOTTING THE RESIDUALS FOR THE BEST FIT TO THE DATA
subplot(4,1,4); hold on;
bar(X, R, 'facecolor', [0 0 0]);
gca_props(); grid on;
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
xlabel('$$ \bf  k_z (A^{-1}) $$', 'Interpreter', 'latex');

end