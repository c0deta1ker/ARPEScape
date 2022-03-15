function fig = kz2model_view_init(xdat, ydat, cTYPE, iparams)
% fig = kz2model_view_init(xdat, ydat, cTYPE, iparams)
%   This function is used to plot the initial curve fitting model prior to
%   using the 'kz2model_solver()' algorithm. This is used as an informative 
%   plot that allows you to view and create a better initial guess of the
%   model prior to running the fitting algorithm.
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
X   = xdat;                 % Original domain
D   = ydat;                 % Original data
M0  = PESCurve(xdat, cTYPE, iparams{1}(1), iparams{1}(2), 1./iparams{1}(3), iparams{1}(5), iparams{1}(6), iparams{1}(7), 0.0, iparams{1}(8));                   % Model data (MFP only)
M1  = PESCurve(xdat, cTYPE, mean([iparams{1}(1),iparams{1}(1)+iparams{1}(6)]), iparams{1}(2), 1./iparams{1}(4), iparams{1}(5), 0.0, 0.0, 0.0, iparams{1}(8));                   % Model data (KZ only)
M   = PESCurve(xdat, cTYPE, iparams{1}(1), iparams{1}(2), 1./iparams{1}(3)+1./iparams{1}(4), iparams{1}(5), iparams{1}(6), iparams{1}(7), 0.0, iparams{1}(8));	% Model data (MFP+KZ)
B   = iparams{1}(9);        % Background
DB	= D - B;                % Data after being background subtracted
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
area(X, M, 'FaceColor', pp.col.fit{2}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
plot(X, M, 'k-', 'linewidth', 0.25);
area(X, M0, 'FaceColor', pp.col.fit{3}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
plot(X, M0, 'k-', 'linewidth', 0.25);
area(X, M1, 'FaceColor', pp.col.fit{4}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
plot(X, M1, 'k-', 'linewidth', 0.25);
% -- Plotting the curve component energy locations
line([iparams{1}(1), iparams{1}(1)], [0, iparams{1}(2)], 'Color', [0 0 0], 'LineWidth', 2, 'Linestyle', '-');
line(iparams{1}(1)+[iparams{1}(6), iparams{1}(6)], [0, iparams{1}(2)*iparams{1}(7)], 'Color', [0 0 0], 'LineWidth', 2, 'Linestyle', '-');
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
subplot(4,1,4); hold on;
bar(X, R, 'facecolor', [0 0 0]);
gca_props(); grid on;
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
xlabel('$$ \bf  k_z (A^{-1}) $$', 'Interpreter', 'latex');

%% - 5 - Plotting the uncertainties to allow for easier optimisation
subplot(4,1,[1,2,3]); hold on;
%% - 5.1 - Plotting the primary peak uncertainties
lbBE  	= iparams{2}(1); ubBE      = iparams{3}(1);
lbINT  	= iparams{2}(2); ubINT     = iparams{3}(2);
x_vals = [lbBE, lbBE, ubBE, ubBE, lbBE];
y_vals = [lbINT, ubINT, ubINT, lbINT, lbINT];
patch(x_vals, y_vals, pp.col.fit{3}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);

end