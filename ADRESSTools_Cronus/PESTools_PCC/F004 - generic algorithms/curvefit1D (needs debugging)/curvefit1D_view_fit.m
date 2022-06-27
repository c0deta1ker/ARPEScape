function fig = curvefit1D_view_fit(fitStr)
% fig = curvefit1D_view_fit(fitStr)
%   This function is used to plot the results of the curve fitting
%   performed by 'fitStr = curvefit1D_solver(...)'. 
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   fitStr:     MATLAB data-structure that contains all fit parameters / variables / information
%
%   OUT:
%   -   fig:    	MATLAB figure object showing the final best fit

%% Initialising the plot parameters
pp  = plot_props();
% - Finding the optimal y-limits
minY    = min([min(fitStr.cYY(:)), min(fitStr.M(:)), min(fitStr.D(:)), min(fitStr.B(:))]);
maxY    = 1.25*max([max(fitStr.cYY(:)), max(fitStr.M(:)), max(fitStr.D(:)), max(fitStr.B(:))]);
ylims   = [minY, maxY];
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'Curve Fit 1D - Final Fit');
fig.Position(3) = 2.5*pp.fig5x4(1); 
fig.Position(4) = 1.0*pp.fig5x4(2);
% -- Plotting the raw data and background
subplot(121); hold on;
h = patch(...
    [fitStr.init_bgrnd{1}(1), fitStr.init_bgrnd{1}(1), fitStr.init_bgrnd{1}(2), fitStr.init_bgrnd{1}(2), fitStr.init_bgrnd{1}(1)],...
    [-1, 1, 1, -1, -1].*1e4, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', [0 0 0]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(fitStr.xdat, fitStr.ydat, 'b-', 'linewidth', 0.5*pp.lwidth);
plot(fitStr.X, fitStr.D, 'b-', 'linewidth', pp.llwidth);
plot(fitStr.X, fitStr.B, 'r-', 'linewidth', pp.llwidth);
plot(fitStr.X, fitStr.DB, 'k-', 'linewidth', pp.llwidth);
gca_props(); grid on;
xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
axis([min(fitStr.xdat(:)), max(fitStr.xdat(:)), ylims]);
legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
title("Background Subtraction", 'interpreter', 'none', 'fontsize', 9);
% -- Plotting the best fit curve components and the initial model
subplot(4,2,[2,4,6]); hold on;
% -- Plotting all of the curve components
for i = 1:fitStr.ncurves
    area(fitStr.XX, fitStr.cYY(:,i), 'FaceColor', pp.col.fit{i}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
    plot(fitStr.XX, fitStr.cYY(:,i), 'k-', 'linewidth', 0.25);
end
% -- Plotting the curve component energy locations
for i = 1:fitStr.ncurves
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    line([fitStr.init_params{1}(i,1), fitStr.init_params{1}(i,1)], [0, max(fitStr.cYY(:,i))],...
        'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '-');
end
% -- Plotting the experimental and fit spectra
plot(fitStr.X, fitStr.DB, 'k-', 'color', pp.col.dat{1}, 'linewidth', 2*pp.llwidth);
plot(fitStr.X, fitStr.M, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Add annotation for the quality of fit
text(0.04, 0.90, "$$ \chi^2 = $$ " + string(fitStr.CHISQ),...
    'interpreter', 'latex', 'fontsize', 14, 'color', 'k', 'Units','normalized');
% - Formatting the figure
gca_props(); grid on;
ylabel('$$ \bf  Y$$', 'Interpreter', 'latex');
axis([min(fitStr.X(:)), max(fitStr.X(:)), ylims]);
title("Final Fit", 'interpreter', 'none', 'fontsize', 9);

% -- Plotting the residuals of the model
subplot(4,2,8); hold on;
bar(fitStr.X, fitStr.R, 'facecolor', [0 0 0]);
gca_props(); grid on;
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
