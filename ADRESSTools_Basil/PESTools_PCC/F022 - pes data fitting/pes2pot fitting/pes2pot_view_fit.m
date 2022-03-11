function fig = pes2pot_view_fit(fitStr)
% fig = pes2pot_view_fit(fitStr)
%   This function is used to plot the initial, model XPS curve PRIOR to
%   curve fitting with 'xps_solver()'. The plot consists of 3 subplots; (1) The
%   background that is determined from the fit; (2) A plot showing all of
%   the fitted curve components, as well as the final model fit and
%   experimental data; (3) A plot of the residuals, showing the quality of
%   the experimental and model fit. This is used as an informative plot
%   that allows you to view and create a better initial guess of the XPS
%   model prior to running the fitting algorithm.
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   pesStr:         MATLAB data-structure that contains the XPS data.
%   -   cTYPE:          1xN vector of the type of curve to use for the nth state.
%   -   iPESCurves:   	3 cells {x0}{lb}{ub} with Nx8 arrays: the n'th peak parameters [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%   -   iPESBgrnd:      3 cells {x0}{lb}{ub} with 1x5 vectors: the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%
%   OUT:
%   -   fig:  	MATLAB figure object with the ARPES data plotted.


%% - 1 - Plotting the model to be used for fitting
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'XPS Curve Fitting');
fig.Position(1) = 100;
fig.Position(2) = 100;
fig.Position(3) = 2.5*pp.fig5x4(1); 
fig.Position(4) = 2.5*pp.fig5x4(2);

%% - 2.1 - PLOTTING THE RAW DATA AND BEST FIT BACKGROUND
subplot(221); hold on;
h = patch(...
    [fitStr.ibgrnd{1}(1), fitStr.ibgrnd{1}(1), fitStr.ibgrnd{1}(2), fitStr.ibgrnd{1}(2), fitStr.ibgrnd{1}(1)],...
    [-1, 1, 1, -1, -1].*1e4, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', [0 0 0]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(fitStr.xdat, fitStr.ydat, 'b-', 'linewidth', 0.5*pp.lwidth);
plot(fitStr.X, fitStr.D, 'b-', 'linewidth', pp.llwidth);
plot(fitStr.X, fitStr.B, 'r-', 'linewidth', pp.llwidth);
plot(fitStr.X, fitStr.DB, 'k-', 'linewidth', pp.llwidth);
gca_props(); grid on;
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
axis([min(fitStr.xdat(:)), max(fitStr.xdat(:)),0, 1.25*max(fitStr.ydat(:))]);
legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
title("Background Subtraction", 'interpreter', 'none', 'fontsize', 9);

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
subplot(223); hold on;
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
img_props(); colormap jet; title('Model Potential : V(z)', 'Interpreter','none');
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