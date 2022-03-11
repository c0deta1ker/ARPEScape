function fig = pes2pot_view_init(pesStr, cTYPE, iparams, bTYPE, ibgrnd, MFP, ZPOT, EPOT)
% fig = pes2pot_view_init(pesStr, cTYPE, iparams, bTYPE, ibgrnd, MFP, ZPOT, EPOT)
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

%% Extracting the original data and interpolating it to a finer grid size
xdat	= linspace(min(pesStr.xdat(:)), max(pesStr.xdat(:)), 1e3);
ydat	= interp1(pesStr.xdat, pesStr.ydat, xdat);

%% - 1 - Extracting all background subtraction data and information
% -- Extracting the background
[X, D, B]  = PESBackground(xdat, ydat,...
    bTYPE, ibgrnd{1}(1), ibgrnd{1}(2), ibgrnd{1}(3), ibgrnd{1}(4), ibgrnd{1}(5), ibgrnd{1}(6));

%% - 2 - Extracting all model data and information
% -- Extracting the total number of states to be fitted
nSTATES = length(cTYPE);
BOFF  	= -1*EPOT(1);
% -- Initialising the vectors to contain XPS data
DB    = D - B;              % Data after being background subtracted
M     = zeros(size(X));   	% Initialising model data
% -- Filing through all curve components and extracting them seperately
cYY = []; potYY = [];
for i = 1:nSTATES
    % --- Extracting curve profiles
   	[cYY(:,i), potYY(:,:,i)] = PESCurve_POT(X, cTYPE(i),...
        iparams{1}(i,1), iparams{1}(i,2), iparams{1}(i,3),...
        iparams{1}(i,4), iparams{1}(i,5), iparams{1}(i,6),...
        iparams{1}(i,7), iparams{1}(i,8), MFP, ZPOT, EPOT);
    M = M + cYY(:,i);
    % --- Storing each curve parameter
    DCL(i)      = iparams{1}(i,1);
    INT(i)      = iparams{1}(i,2);
    FWHM(i)     = iparams{1}(i,3);
    MR(i)       = iparams{1}(i,4);
    LSE(i)      = iparams{1}(i,5);
    LSI(i)      = iparams{1}(i,6);
    LSW(i)      = iparams{1}(i,7);
    ASY(i)      = iparams{1}(i,8);
    % --- Storing additional parameters
    BE_Z0(i)	= DCL(i) + EPOT(1);
    [~,I] = max(cYY(:,i)); BE_MAX(i)= X(I);
end
% - Extracting the total potential profile of all components
potYY_tot         = sum(potYY, 3);
potYY_comp        = squeeze(sum(potYY, 2));
for i = 1:nSTATES
    potYY_comp(:,i) = INT(i) .* (potYY_comp(:,i) ./ max(potYY_comp(:,i)));
end

%% - 3 - Determination of the residuals and chi-squared
R           = M - (D - B);                  % Residuals
CHISQ       = sum(R.^2 ./ abs(M));          % Chi-squared

%% - 4 - Plotting the model to be used for fitting
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'XPS Curve Fitting');
fig.Position(1) = 100;
fig.Position(2) = 100;
fig.Position(3) = 2.5*pp.fig5x4(1); 
fig.Position(4) = 2.5*pp.fig5x4(2);

%% - 4.1 - PLOTTING THE RAW DATA AND BEST FIT BACKGROUND
subplot(221); hold on;
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
subplot(8,2,[2,4,6]); hold on;
% -- Plotting all of the curve components
for i = 1:nSTATES
    area(X, cYY(:,i), 'FaceColor', pp.col.fit{i}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
    plot(X, cYY(:,i), 'k-', 'linewidth', 0.25);
end
% -- Plotting the curve component energy locations
for i = 1:nSTATES
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    % --- Plotting DCL
    line([1, 1]*DCL(i), [0, max(cYY(:,i))], 'Color', pp.col.fit{i}, 'LineWidth', lWidth, 'Linestyle', '-');
    % --- Plotting BE
    line([1, 1]*DCL(i)+EPOT(1), [0, max(cYY(:,i))], 'Color', pp.col.fit{i}, 'LineWidth', lWidth, 'Linestyle', ':');
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
title("Initial Model for Fitting", 'interpreter', 'none', 'fontsize', 9);

%% - 4.3 - PLOTTING THE RESIDUALS FOR THE BEST FIT TO THE DATA
subplot(8,2,8); hold on;
bar(X, R, 'facecolor', [0 0 0]);
gca_props(); grid on;
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');

%% - 5 - Plotting the uncertainties to allow for easier optimisation
subplot(8,2,[2,4,6]); hold on;
for i = 1:nSTATES
    %% - 5.1 - Plotting the primary peak uncertainties
    lbBE  	= iparams{2}(i,1); ubBE      = iparams{3}(i,1);
    lbINT  	= iparams{2}(i,2); ubINT     = iparams{3}(i,2);
    x_vals = [lbBE, lbBE, ubBE, ubBE, lbBE];
    y_vals = [lbINT, ubINT, ubINT, lbINT, lbINT];
    patch(x_vals, y_vals, pp.col.fit{i}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);
    %% - 5.2 - Adding text to identify each peak index
    text(max(x_vals), max(y_vals), string(i), 'interpreter', 'latex', 'fontsize', 14, 'color', 'k');
    %% - 5.3 - Plotting the spin-orbit split (SOS) peak uncertainties
    DCL_mu  = DCL(i);
    lbLSE  	= iparams{2}(i,5);   ubLSE	= iparams{3}(i,5);
    lbLSI  	= iparams{2}(i,6);   ubLSI	= iparams{3}(i,6);
    x_vals  = DCL_mu + [lbLSE, lbLSE, ubLSE, ubLSE, lbLSE];
    y_vals  = [lbINT*lbLSI, ubINT*ubLSI, ubINT*ubLSI, lbINT*lbLSI, lbINT*lbLSI];
    patch(x_vals, y_vals, pp.col.fit{i}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);
end

%% - 6 - PLOTTING THE POTENTIAL WELL, MODULATED BY THE MFP
subplot(223); hold on;
% --- Image of the curve series shifted by potential and scaled by MFP
ImData(ZPOT, X, potYY_tot);
% -- Plotting the curve component energy locations
for i = 1:nSTATES
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    % --- Plotting DCL
    line([0, max(ZPOT(:))], [DCL(i), DCL(i)], 'Color', pp.col.fit{i}, 'LineWidth', lWidth, 'Linestyle', '-');
    % --- Plotting BE
    line([0, max(ZPOT(:))], [DCL(i), DCL(i)]+ EPOT(1), 'Color', pp.col.fit{i}, 'LineWidth', lWidth, 'Linestyle', ':');
end
% --- Plotting the potential energy curves
plot(ZPOT, EPOT+DCL(1), 'r-', 'linewidth', 2);
% --- Formatting the figure
img_props(); colormap jet; title('Model Potential : V(z)', 'Interpreter','none');
axis square;
xlabel('$$ \bf  z (nm) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
axis([0, max(ZPOT(:)), min(X(:)), max(X(:))]);
% --- Adding text to quantify the band offset
BOFF        = round(BOFF, 3);
text(0.04, 0.92, "$$ \phi_{z=0} = $$ " + string(BOFF), 'interpreter', 'latex', 'fontsize', 14, 'color', 'w', 'Units','normalized');

%% - 7 - PLOTTING ALL THE INDIVIDUAL CURVES, SHIFTED BY THE POTENTIAL FOR N=1
subplot(224); hold on;
% --- Plotting individual curves
cols = flipud(jet(length(ZPOT)+2));
for i = 1:length(ZPOT)
    plot(X, potYY(:,i,1), 'k-', 'color', cols(i,:), 'linewidth', 0.5);
end
% --- Plotting the final curve
plot(X, potYY_comp(:,1), 'k-', 'linewidth', 2);
% --- Formatting the figure
gca_props(); title('Integrated Curves : n = 1', 'Interpreter','none');
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
axis([min(X(:)), max(X(:)), 0, 1.1*max(M(:))]);
% -- Plotting the curve component energy locations
lWidth = pp.llwidth;
% --- Plotting DCL
line([1, 1]*DCL(1), [0, max(cYY(:,1))], 'Color', pp.col.fit{1}, 'LineWidth', lWidth, 'Linestyle', '-');
% --- Plotting BE
line([1, 1]*DCL(1)+EPOT(1), [0, max(cYY(:,1))], 'Color', pp.col.fit{1}, 'LineWidth', lWidth, 'Linestyle', ':');
% --- Plotting BE_MAX
[~,I] = max(potYY_comp(:,1)); BE_MAX = X(I);
line([BE_MAX, BE_MAX], [0, max(cYY(:,1))], 'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '--');
% --- Adding text to quantify all the shifts
DCL         = round(DCL(1), 3);
BE_Z0       = round(BE_Z0(1), 3);
BE_MAX      = round(BE_MAX(1), 3);
text(0.04, 0.92, "$$ \Delta_{CL} = $$ " + string(DCL), 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');
text(0.04, 0.86, "$$ BE_{z=0} = $$ " + string(BE_Z0), 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');
text(0.04, 0.80, "$$ BE_{MAX} = $$ " + string(BE_MAX), 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');

end