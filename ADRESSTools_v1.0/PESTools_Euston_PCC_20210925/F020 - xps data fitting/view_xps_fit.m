function fig = view_xps_fit(fitStr, show_bgrnd)
% fig = view_xps_fit(fitStr, show_bgrnd)
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
if nargin < 2; show_bgrnd = 1; end
if isempty(show_bgrnd); show_bgrnd = 1; end
% -- Initialising the plot parameters
pp  = plot_props();

%% - 1 - FIGURE WITH BACKGROUND AND XPS FITS
if show_bgrnd == 1
    %% - 1.1 - INITIALISING THE FIGURE
    fig = figure(); set(fig, 'Name', 'XPS Curve Fitting');
    fig.Position(3) = 2.5*pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    %% - 1.2 - PLOTTING THE RAW DATA AND BEST FIT BACKGROUND
    subplot(121); hold on;
    h = patch(...
        [min(fitStr.X(:)), min(fitStr.X(:)), max(fitStr.X(:)), max(fitStr.X(:)), min(fitStr.X(:))],...
    [-1, 1, 1, -1, -1].*1e4, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', [0 0 0]);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot(fitStr.xdat, fitStr.int, 'b-', 'linewidth', 0.5*pp.lwidth);
    plot(fitStr.X, fitStr.D, 'b-', 'linewidth', pp.llwidth);
    plot(fitStr.X, fitStr.B, 'r-', 'linewidth', pp.llwidth);
    plot(fitStr.X, fitStr.DB, 'k-', 'linewidth', pp.llwidth);
    gca_props(); grid on;
    xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
    if min(fitStr.DB(:)) > 0; min_y_val = 0; 
    else; min_y_val = min(fitStr.DB(:)); 
    end
    axis([min(fitStr.xdat(:)), max(fitStr.xdat(:)), min_y_val, 1.25*max(fitStr.int(:))]);
    legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
    title("Background Subtraction", 'interpreter', 'none', 'fontsize', 9);
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
    
%% - 2 - FIGURE WITH XPS FIT ONLY
else
    %% - 2.1 - INITIALISING THE FIGURE
    fig = figure(); set(fig, 'Name', 'XPS Curve Fitting');
    fig.Position(3) = pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2); 
    %% - 2.2 - PLOTTING THE RESIDUALS FOR THE BEST FIT TO THE DATA
    subplot(3,1,3); hold on;
    bar(fitStr.X, fitStr.R, 'facecolor', [0 0 0]);
    gca_props(); grid on;
    xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
    %% - 2.3 - PLOTTING THE BEST FIT CURVE COMPONENTS AND FINAL MODEL FIT
    subplot(3,1,[1,2]); hold on;
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
    text(0.04, 0.92, "$$ \chi^2 = $$ " + string(fitStr.CHISQ),...
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
end
    
