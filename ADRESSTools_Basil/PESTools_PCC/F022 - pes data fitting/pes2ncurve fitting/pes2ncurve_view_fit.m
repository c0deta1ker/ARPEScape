function fig = pes2ncurve_view_fit(fitStr)
% fig = pes2ncurve_view_fit(fitStr, show_bgrnd)
%   This function is used to plot the results of the XPS curve fitting
%   performed by 'xps_solver()'. The plot consists of 3 subplots; (1) The
%   background that is determined from the fit; (2) A plot showing all of
%   the fitted curve components, as well as the final model fit and
%   experimental data; (3) A plot of the residuals, showing the quality of
%   the experimental and model fit. If you only want to save the fitted XPS
%   spectrum and the residuals (without the background subtraction
%   discussed in (1), you can set the 'show_bgrnd' variable to 0; its
%   default value is 1.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   fitStr:	data    structure of the XPS data after fitting.
%   -   show_bgrnd:     either 1 or 0; 1 to show background and XPS fits, 0 to only show XPS fits.
%
%   OUT:
%   -   fig:            MATLAB figure object with the ARPES data plotted.

%% Default parameters
% -- Initialising the plot parameters
pp  = plot_props();

%% - 1 - FIGURE WITH BACKGROUND AND XPS FITS
%% - 1.1 - INITIALISING THE FIGURE
fig = figure(); set(fig, 'Name', 'PES Curve Fitting');
fig.Position(3) = 2.5*pp.fig5x4(1); 
fig.Position(4) = pp.fig5x4(2);
%% - 1.2 - PLOTTING THE RAW DATA AND BEST FIT BACKGROUND
subplot(121); hold on;
% -- Plotting the ROI, LHS and RHS analysis windows
ROI     = [fitStr.X(1), fitStr.X(end)]; 
LHS     = ROI(1) + 0.05.*[-1,1].*range(fitStr.xdat(:));
RHS     = ROI(2) + 0.05.*[-1,1].*range(fitStr.xdat(:));
hLHS    = patch([LHS(1), LHS(1), LHS(2), LHS(2), LHS(1)], [-1, 1, 1, -1, -1].*1e6, [0.6 0.6 0.2], 'facealpha', 0.4, 'edgecolor', 'none');
hLHS.Annotation.LegendInformation.IconDisplayStyle = 'off';
hRHS    = patch([RHS(1), RHS(1), RHS(2), RHS(2), RHS(1)], [-1, 1, 1, -1, -1].*1e6, [0.6 0.6 0.2], 'facealpha', 0.4, 'edgecolor', 'none');
hRHS.Annotation.LegendInformation.IconDisplayStyle = 'off';
hMID    = patch([ROI(1), ROI(1), ROI(2), ROI(2), ROI(1)], [-1, 1, 1, -1, -1].*1e6, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', 'none');
hMID.Annotation.LegendInformation.IconDisplayStyle = 'off';
% -- Plotting a vertical line to show the ROI
a = line([ROI(1) ROI(1)], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 1, 'Linestyle', '-');
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
b = line([ROI(2) ROI(2)], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 1, 'Linestyle', '-');
b.Annotation.LegendInformation.IconDisplayStyle = 'off';
% -- Plotting the 1D data
plot(fitStr.xdat, fitStr.ydat, 'b-', 'linewidth', 0.5);
plot(fitStr.X, fitStr.D, 'b-', 'linewidth', 2);
plot(fitStr.X, fitStr.B, 'r-', 'linewidth', 2);
plot(fitStr.X, fitStr.DB, 'k-', 'linewidth', 2);
gca_props(); title('Background Subtraction', 'interpreter', 'none', 'fontsize', 9); 
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity$$', 'Interpreter', 'latex');
legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
% -- Determining the best limits for the plot
axLim_y = [fitStr.ydat; fitStr.DB];
axis([min(fitStr.xdat(:)), max(fitStr.xdat(:)), min(axLim_y(:)), 1.1*max(axLim_y(:))]);
%% - 1.3 - PLOTTING THE RESIDUALS FOR THE BEST FIT TO THE DATA
subplot(4,2,8); hold on;
bar(fitStr.X, fitStr.R, 'facecolor', [0 0 0]);
gca_props(); grid on;
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
%% - 1.4 - PLOTTING THE BEST FIT CURVE COMPONENTS AND FINAL MODEL FIT
subplot(4,2,[2,4,6]); hold on;
% -- Plotting all of the curve components
for i = 1:fitStr.nSTATES
    area(fitStr.XX, fitStr.cYY(:,i), 'FaceColor', pp.col.fit{i}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
    plot(fitStr.XX, fitStr.cYY(:,i), 'k-', 'linewidth', 0.25);
end
% -- Plotting the curve component energy locations
for i = 1:fitStr.nSTATES
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    line([fitStr.cPARAMS(i,1), fitStr.cPARAMS(i,1)], [0, max(fitStr.cYY(:,i))],...
        'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '-');
end
% -- Plotting the experimental and fit spectra
plot(fitStr.X, fitStr.DB, 'k-', 'color', pp.col.dat{1}, 'linewidth', 2*pp.llwidth);
plot(fitStr.XX, fitStr.YY, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Add annotation for the quality of fit
text(0.04, 0.90, "$$ \chi^2 = $$ " + string(fitStr.CHISQ),...
    'interpreter', 'latex', 'fontsize', 14, 'color', 'k', 'Units','normalized');
% - Formatting the figure
gca_props(); grid on;
ylabel('$$ \bf  Intensity$$', 'Interpreter', 'latex');
if min(fitStr.DB(:)) < fitStr.YY; min_y_val = min(fitStr.DB(:));
else; min_y_val = min(fitStr.YY(:)); 
end
if min_y_val > 0; min_y_val = 0; end
axis([min(fitStr.X(:)), max(fitStr.X(:)), min_y_val, 1.10*max(fitStr.DB(:))]);
title("Best Fit PES Curves", 'interpreter', 'none', 'fontsize', 9);
