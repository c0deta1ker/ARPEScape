function fig = pes2boff_view_fit(fitStr)
% fig = pes2boff_view_fit(fitStr)
%   This function is used to plot the final PES curve model AFTER  curve 
%   fitting with the 'pes2boff_solver()'.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   fitStr: data structure that contains the PES fitting output of the corresponding solver.
%
%   OUT:
%   -   fig:    MATLAB figure object that summarises the fit.

%% - 1 - Plotting the model to be used for fitting
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'PES Curve Fitting');
fig.Position(1) = 100;
fig.Position(2) = 100;
fig.Position(3) = 2.5*pp.fig5x4(1); 
fig.Position(4) = 2.5*pp.fig5x4(2);

%% - 2.1 - PLOTTING THE RAW DATA AND BEST FIT BACKGROUND
subplot(221); hold on;
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

%% - 2.2 - PLOTTING THE BEST FIT CURVE COMPONENTS AND FINAL MODEL FIT
subplot(8,2,[2,4,6]); hold on;
% -- Plotting all of the curve components
for i = 1:fitStr.nSTATES
    area(fitStr.XX, fitStr.cYY(:,i), 'FaceColor', pp.col.fit{i}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
    plot(fitStr.XX, fitStr.cYY(:,i), 'k-', 'linewidth', 0.25);
end
% -- Plotting the curve component energy locations
for i = 1:fitStr.nSTATES
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    % --- Plotting DCL
    line([1,1]*fitStr.DCL(i), [0, max(fitStr.cYY(:,i))], 'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '-');
    % --- Plotting BE_Z0
    line([1,1]*fitStr.BE_Z0(i), [0, max(fitStr.cYY(:,i))], 'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', ':');
end
% -- Plotting the experimental and fit spectra
plot(fitStr.X, fitStr.DB, 'k-', 'color', pp.col.dat{1}, 'linewidth', 2*pp.llwidth);
plot(fitStr.X, fitStr.M, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Add annotation for the quality of fit
text(0.04, 0.90, "$$ \chi^2 = $$ " + string(fitStr.CHISQ),...
    'interpreter', 'latex', 'fontsize', 14, 'color', 'k', 'Units','normalized');
% - Formatting the figure
gca_props(); grid on;
ylabel('$$ \bf  Intensity$$', 'Interpreter', 'latex');
axis([min(fitStr.X(:)), max(fitStr.X(:)), min(fitStr.DB(:)), 1.10*max(fitStr.DB(:))]);
title("Initial Model for Fitting", 'interpreter', 'none', 'fontsize', 9);

%% - 2.3 - PLOTTING THE RESIDUALS FOR THE BEST FIT TO THE DATA
subplot(8,2,8); hold on;
bar(fitStr.X, fitStr.R, 'facecolor', [0 0 0]);
gca_props(); grid on;
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');

%% - 3 - PLOTTING THE POTENTIAL WELL, MODULATED BY THE MFP
subplot(8,2,[9,11,13,15]); hold on;
% --- Image of the curve series shifted by potential and scaled by MFP
ImData(fitStr.ZPOT, fitStr.XX, fitStr.potYY_tot);
% -- Plotting the curve component energy locations
for i = 1:fitStr.nSTATES
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    % --- Plotting DCL
    line([0, max(fitStr.ZPOT(:))], [1,1]*fitStr.DCL(i), 'Color', pp.col.fit{i}, 'LineWidth', lWidth, 'Linestyle', '-');
    % --- Plotting BE_Z0
    line([0, max(fitStr.ZPOT(:))], [1,1]*fitStr.BE_Z0(i), 'Color', pp.col.fit{i}, 'LineWidth', lWidth, 'Linestyle', ':');
end
% --- Plotting the potential energy curves
plot(fitStr.ZPOT, fitStr.EPOT+fitStr.DCL(1), 'r-', 'linewidth', 2);
% --- Formatting the figure
img_props(); 
subplot(8,2,[9,11,13,15]); hold on;
colormap jet; title('Model Potential : V(z)', 'Interpreter','none');
axis square;
xlabel('$$ \bf  z (nm) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
axis([0, max(fitStr.ZPOT(:)), min(fitStr.XX(:)), max(fitStr.XX(:))]);
% --- Adding text to quantify the band offset
BOFF        = round(fitStr.BOFF, 3);
text(0.04, 0.92, "$$ \phi_{z=0} = $$ " + string(BOFF), 'interpreter', 'latex', 'fontsize', 14, 'color', 'w', 'Units','normalized');

%% - 4 - PLOTTING ALL THE INDIVIDUAL CURVES, SHIFTED BY THE POTENTIAL FOR N=1
subplot(224); hold on;
% --- Plotting individual curves
cols = flipud(jet(length(fitStr.ZPOT)+2));
for i = 1:length(fitStr.ZPOT)
    plot(fitStr.XX, fitStr.potYY(:,i,1), 'k-', 'color', cols(i,:), 'linewidth', 0.5);
end
% --- Plotting the final curve
plot(fitStr.XX, fitStr.potYY_comp(:,1), 'k-', 'linewidth', 2);
% --- Formatting the figure
gca_props(); title('Integrated Curves : n = 1', 'Interpreter','none');
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
axis([min(fitStr.XX(:)), max(fitStr.XX(:)), 0, 1.1*max(fitStr.M(:))]);
% -- Plotting the curve component energy locations
lWidth = pp.llwidth;
% --- Plotting DCL
line([1,1]*fitStr.DCL(1), [0, max(fitStr.cYY(:,1))], 'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '-');
% --- Plotting BE_Z0
line([1,1]*fitStr.BE_Z0(1), [0, max(fitStr.cYY(:,1))], 'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', ':');
% --- Plotting BE_MAX
[~,I] = max(fitStr.potYY_comp(:,1)); BE_MAX = fitStr.XX(I);
line([BE_MAX, BE_MAX], [0, max(fitStr.cYY(:,1))], 'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '--');
% --- Adding text to quantify all the shifts
DCL         = round(fitStr.DCL(1), 3);
BE_Z0       = round(fitStr.BE_Z0(1), 3);
BE_MAX      = round(fitStr.BE_MAX(1), 3);
text(0.04, 0.92, "$$ \Delta_{CL} = $$ " + string(DCL), 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');
text(0.04, 0.86, "$$ BE_{z=0} = $$ " + string(BE_Z0), 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');
text(0.04, 0.80, "$$ BE_{MAX} = $$ " + string(BE_MAX), 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');

end