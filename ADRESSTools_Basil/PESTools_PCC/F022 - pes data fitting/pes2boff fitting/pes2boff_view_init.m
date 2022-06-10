function fig = pes2boff_view_init(pesStr, modelStr, cTYPE, iparams, bTYPE, ibgrnd)
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

%% - 1 - Extracting all background subtraction data and information
% -- Extracting the background
[X, D, B]  = PESBackground(pesStr.xdat, pesStr.ydat,...
    bTYPE, ibgrnd{1}(1), ibgrnd{1}(2), ibgrnd{1}(3), ibgrnd{4});

%% - 2 - Extracting all of the input parameters
% -- Extracting the original XPS data
% xdat	= linspace(min(pesStr.xdat(:)), max(pesStr.xdat(:)), 1e3);
% ydat	= interp1(pesStr.xdat, pesStr.ydat, xdat);
xdat   	= pesStr.xdat;
ydat 	= pesStr.ydat;
% -- Extracting the total number of states to be fitted
nSTATES = length(cTYPE);
% -- Extracting all of the fit parameters
% -- Singular variables
MFP     = iparams{1}(1);
BOFF	= iparams{1}(2);
lbMFP  	= iparams{2}(1);  ubMFP     = iparams{3}(1);
lbBOFF	= iparams{2}(2);  ubBOFF    = iparams{3}(2);
% -- Curve parameter variables
for i = 1:nSTATES
    % --- Storing each curve parameter
    DCL(i)      = iparams{1}(i+2+0*nSTATES);
    INT(i)      = iparams{1}(i+2+1*nSTATES);
    FWHM(i)     = iparams{1}(i+2+2*nSTATES);
    MR(i)       = iparams{1}(i+2+3*nSTATES);
    LSE(i)      = iparams{1}(i+2+4*nSTATES);
    LSI(i)      = iparams{1}(i+2+5*nSTATES);
    LSW(i)      = iparams{1}(i+2+6*nSTATES);
    ASY(i)      = iparams{1}(i+2+7*nSTATES);
    % --- Storing the constraints of each curve parameter
    lbDCL(i)  	= iparams{2}(i+2+0*nSTATES);  ubDCL(i)       = iparams{3}(i+2+0*nSTATES);
    lbINT(i)  	= iparams{2}(i+2+1*nSTATES);  ubINT(i)       = iparams{3}(i+2+1*nSTATES);
    lbFWHM(i)	= iparams{2}(i+2+2*nSTATES);  ubFWHM(i)      = iparams{3}(i+2+2*nSTATES);
    lbMR(i)  	= iparams{2}(i+2+3*nSTATES);  ubMR(i)        = iparams{3}(i+2+3*nSTATES);
    lbLSE(i)  	= iparams{2}(i+2+4*nSTATES);  ubLSE(i)       = iparams{3}(i+2+4*nSTATES);
    lbLSI(i)  	= iparams{2}(i+2+5*nSTATES);  ubLSI(i)       = iparams{3}(i+2+5*nSTATES);
    lbLSW(i)  	= iparams{2}(i+2+6*nSTATES);  ubLSW(i)       = iparams{3}(i+2+6*nSTATES);
    lbASY(i)  	= iparams{2}(i+2+7*nSTATES);  ubASY(i)       = iparams{3}(i+2+7*nSTATES);
    % --- Storing additional parameters
    BE_Z0(i)	= DCL(i) - BOFF;
end
% -- Extracting the model potential well based on the model and BOFF
imodelStr     	= mstheory_interp(modelStr, BOFF, MFP);
ZPOT            = linspace(min(imodelStr.ZPOT(:)), max(imodelStr.ZPOT(:)), 1e3);
EPOT            = interp1(imodelStr.ZPOT, imodelStr.EPOT, ZPOT);

%% - 3 - Extracting all model data and information
% -- Initialising the vectors to contain XPS data
DB    = D - B;              % Data after being background subtracted
M     = zeros(size(X));   	% Initialising model data
% -- Filing through all curve components and extracting them seperately
cYY = []; potYY = [];
for i = 1:nSTATES
    % --- Extracting curve profiles
   	[cYY(:,i), potYY(:,:,i)] = PESCurve_POT(X, cTYPE(i),...
        DCL(i), INT(i), FWHM(i), MR(i), LSE(i), LSI(i), LSW(i), ASY(i),...
        MFP, ZPOT, EPOT);
    M = M + cYY(:,i);
    % --- Storing additional parameters
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
CHISQ       = sum(R.^2 ./ abs(M + B));      % Chi-squared

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
% -- Plotting the ROI, LHS and RHS analysis windows
ROI     = [X(1), X(end)]; 
LHS     = ROI(1) + 0.05.*[-1,1].*range(xdat(:));
RHS     = ROI(2) + 0.05.*[-1,1].*range(xdat(:));
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
plot(xdat, ydat, 'b-', 'linewidth', 0.5);
plot(X, D, 'b-', 'linewidth', 2);
plot(X, B, 'r-', 'linewidth', 2);
plot(X, DB, 'k-', 'linewidth', 2);
gca_props(); title('Background Subtraction', 'interpreter', 'none', 'fontsize', 9); 
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity$$', 'Interpreter', 'latex');
legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
% -- Determining the best limits for the plot
axLim_y = [ydat; DB];
axis([min(xdat(:)), max(xdat(:)), min(axLim_y(:)), 1.1*max(axLim_y(:))]);

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
    x_vals = [lbDCL(i), lbDCL(i), ubDCL(i), ubDCL(i), lbDCL(i)];
    y_vals = [lbINT(i), ubINT(i), ubINT(i), lbINT(i), lbINT(i)];
    patch(x_vals, y_vals, pp.col.fit{i}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);
    %% - 5.2 - Adding text to identify each peak index
    text(max(x_vals), max(y_vals), string(i), 'interpreter', 'latex', 'fontsize', 14, 'color', 'k');
    %% - 5.3 - Plotting the spin-orbit split (SOS) peak uncertainties
    x_vals  = DCL(i) + [lbLSE(i), lbLSE(i), ubLSE(i), ubLSE(i), lbLSE(i)];
    y_vals  = [lbINT(i)*lbLSI(i), ubINT(i)*ubLSI(i), ubINT(i)*ubLSI(i), lbINT(i)*lbLSI(i), lbINT(i)*lbLSI(i)];
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
text(0.04, 0.92, "$$ \phi_{z=0} = $$ " + string(BOFF) + " eV", 'interpreter', 'latex', 'fontsize', 14, 'color', 'w', 'Units','normalized');

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
text(0.04, 0.92, "$$ \Delta_{CL} = $$ " + string(DCL) + " eV", 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');
text(0.04, 0.86, "$$ BE_{z=0} = $$ " + string(BE_Z0) + " eV", 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');
text(0.04, 0.80, "$$ BE_{MAX} = $$ " + string(BE_MAX) + " eV", 'interpreter', 'latex', 'fontsize', 13, 'color', 'k', 'Units','normalized');

end