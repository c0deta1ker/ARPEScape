% --- General function to view fitted ARPES data using arpes2model_solver()
function fig = view_arpes2model_fit(fitStr)
% fig = view_arpes2model_fit(fitStr)
%   This function is used to plot the results of the ARPES curve fitting
%   performed by 'arpes2model_solver()'. The plot consists of 3 subplots; (1) The
%   initial model data to be fitted based on the initial parmaeters; (2) A 
%   plot showing all of the experimental ARPES data to be fitted; (3) A 
%   plot summarising the residuals, showing the quality of the experimental
%   and model fit.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   fitStr:  	MATLAB data-structure that is output from 'arpes2model_solver()'.
%
%   OUT:
%   -   fig:    	figure output

%% Initialising input variables
% -- Defining the values of the EDC cut figures
win     = 0.01;
val1    = 0; 
val2    = 0.5*min(fitStr.kx(:));
val3    = 0.5*max(fitStr.kx(:));

%% - 1 - PLOTTING THE FITTED DATA
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'arpes2model fits');
fig.Position(3) = 4.*pp.fig5x4(1);
fig.Position(4) = 2.*pp.fig5x4(2);
%% - 1.1 - Plotting the MODEL data
axs(1) = subplot(3,7,[1,2,8,9]); hold on;
ImData(fitStr.kx, fitStr.eb, fitStr.MB); 
img_props(); 
xlabel(' $$ k_{//} [A^{-1}] $$', 'interpreter', 'latex');
ylabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
axis([fitStr.kx_lims, fitStr.eb_lims]); axis square;
caxis([min(fitStr.MB(:)), max(fitStr.MB(:))]);
text(0.04, 0.94, "$$ \Phi = $$" + string(round(fitStr.BOFF,3)) + "$$ eV $$",...
    'interpreter', 'latex', 'fontsize', 16, 'color', 'w', 'Units','normalized');
title('Best fit model data');

%% - 4.2 - Plotting the ARPES data
axs(2) = subplot(3,7,[3,4,10,11]); hold on;
ImData(fitStr.kx, fitStr.eb, fitStr.D); img_props(); 
xlabel(' $$ k_{//} [A^{-1}] $$', 'interpreter', 'latex');
ylabel('', 'interpreter', 'latex');
axis([fitStr.kx_lims, fitStr.eb_lims]); axis square;
title('ARPES data that is fitted');
caxis([min(fitStr.D(:)), max(fitStr.D(:))]);
% -- Plotting the EDC cut lines
line([1 1].*val1, [-1e5, 1e5], 'Color', pp.col.fit{1}, 'LineWidth', pp.llwidth, 'Linestyle', '-');
line([1 1].*val2, [-1e5, 1e5], 'Color', pp.col.fit{2}, 'LineWidth', pp.llwidth, 'Linestyle', '-');
line([1 1].*val3, [-1e5, 1e5], 'Color', pp.col.fit{3}, 'LineWidth', pp.llwidth, 'Linestyle', '-');

%% - 4.6.1 - LHS and RHS EDCs
axs(7) = subplot(3,7,[17,18]); hold on;
% -- Extracting the LHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.D, 'edc', val2 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.MB, 'edc', val2 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{2}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Extracting the RHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.D, 'edc', val3 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.MB, 'edc', val3 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{3}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
gca_props(); grid on;
xlabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');

%% - 4.3 - Plotting the 2D residuals
axs(3) = subplot(3,7,[5,6,12,13]); hold on;
ImData(fitStr.kx, fitStr.eb, fitStr.CHISQ2D); img_props(); 
axis([fitStr.kx_lims, fitStr.eb_lims]);
title('Chi-Squared');
% -- Add annotation for the quality of fit
text(0.04, 0.94, "$$ \chi^2 = $$ " + string(fitStr.CHISQ),...
    'interpreter', 'latex', 'fontsize', 16, 'color', 'w', 'Units','normalized');
% -- Make the colormap bipolar (negative values blue, positive is red)
colormap(axs(3),'jet');
% -- Make the colormap bipolar (negative values blue, positive is red)
ylabel('', 'interpreter', 'latex');
xlabel('');

%% - 4.4 - Plotting the 1D residuals along x
axs(4) = subplot(3,7,[19,20]); hold on;
bar(fitStr.kx, fitStr.CHISQx);
gca_props(); grid on;
xlabel(' $$ k_{//} [A^{-1}] $$', 'interpreter', 'latex');
ylabel('$$ \bf  \chi^2(x) $$', 'Interpreter', 'latex');
xlim(fitStr.kx_lims);
ax = gca; ax.YAxisLocation = 'right';

%% - 4.5 - Plotting the 1D residuals along y
axs(5) = subplot(3,7,[7,14]); hold on;
barh(fitStr.eb, fitStr.CHISQy);
gca_props(); grid on;
ylabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
xlabel('$$ \bf  \chi^2(y) $$', 'Interpreter', 'latex');
ylim(fitStr.eb_lims);
ax = gca; ax.YAxisLocation = 'right';

%% - 4.6 - Plotting best fit EDCs
%% - 4.6.1 - Gamma EDCs
axs(6) = subplot(3,7,[15,16]); hold on;
% -- Extracting the EDCs
[ARPES_XCut, ARPES_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.D, 'edc', val1 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.MB, 'edc', val1 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{1}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
gca_props(); grid on;
ylabel('$$ Intensity $$', 'Interpreter', 'latex');
xlabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
%% - 4.6.1 - LHS and RHS EDCs
axs(7) = subplot(3,7,[17,18]); hold on;
% -- Extracting the LHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.D, 'edc', val2 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.MB, 'edc', val2 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{2}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Extracting the RHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.D, 'edc', val3 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(fitStr.kx, fitStr.eb, fitStr.MB, 'edc', val3 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k:', 'color', pp.col.fit{3}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k:', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
gca_props(); grid on;
xlabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');

end