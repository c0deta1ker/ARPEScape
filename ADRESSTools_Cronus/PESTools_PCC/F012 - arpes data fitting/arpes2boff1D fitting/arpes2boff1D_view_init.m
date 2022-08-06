function fig = arpes2boff1D_view_init(arpesStr, modelFunc, cTYPE, iparams, bTYPE, ibgrnd)
% fig = arpes2boff1D_view_init(arpesStr, modelFunc, cTYPE, iparams, bTYPE, ibgrnd)
%   This function is used to plot the initial curve fitting model prior to
%   using the 'arpes2boff1D_solver()' algorithm. This is used as an informative 
%   plot that allows you to view and create a better initial guess of the
%   model prior to running the fitting algorithm.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   arpesStr:       MATLAB data-structure that contains the initial ARPES data.
%   -   modelFunc:      1xN cell of functions that gives the band offset vs subband energy
%   -   cTYPE:          1xN vector of the type of curve to use for fitting. Default: "sGLA" ("pGLA", "DS")
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x6 array: the model fit parameters [FDEF,FDT,FDW,BOFF,INT,FWHM]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Poly", "Shir", "LinShir", "StepFDDGpL", "StepFDDGsL")
%   -   ibgrnd:         4 cells {x0}{lb}{ub}{args} with 1x3 vectors of the background parameters: x0=lb=ub=[LHS,RHS,BGR], or argument of the background type args = []
%
%   OUT:
%   -   fig:  	MATLAB figure object with the ARPES data plotted.

%% Default parameters
% - Default input parameters to use
if nargin < 6
    ibgrnd{1} = [...
        mean(arpesStr.xdat(:)) - abs(0.25*range(arpesStr.xdat(:))),...
        mean(arpesStr.xdat(:)) + abs(0.25*range(arpesStr.xdat(:))),...
        0];
    ibgrnd{2} = [0, 0, -abs(0.05*range(arpesStr.ydat(:)))];
    ibgrnd{3} = [0, 0, +abs(0.05*range(arpesStr.ydat(:)))];
    ibgrnd{4} = {1};
end
if nargin < 5; bTYPE = "Poly"; end
if isempty(ibgrnd)
    ibgrnd{1} = [...
        mean(arpesStr.xdat(:)) - abs(0.25*range(arpesStr.xdat(:))),...
        mean(arpesStr.xdat(:)) + abs(0.25*range(arpesStr.xdat(:))),...
        0];
    ibgrnd{2} = [0, 0, -abs(0.05*range(arpesStr.ydat(:)))];
    ibgrnd{3} = [0, 0, +abs(0.05*range(arpesStr.ydat(:)))];
    ibgrnd{4} = {1};
end
if isempty(bTYPE); bTYPE = "Poly"; end
% -- Consistency check and finding the total number of curves
if size(iparams{1}, 1) ~= size(iparams{2}, 1) || size(iparams{2}, 1) ~= size(iparams{3}, 1)
    error('The input parameter cell array is not a consistent size - check iparams input!');
end

%% Initialising variables
% -- Extracting the total number of states to be fitted
nSTATES = length(cTYPE);
% -- Extracting the input variables for the fit
FDEF    = iparams{1}(1);        % scalar of the FDD Fermi-Level position.
FDT     = iparams{1}(2);        % scalar of the FDD temperature.
FDW     = iparams{1}(3);        % scalar of the FDD Gaussian width after convolution.
BOFF    = iparams{1}(4);        % scalar of the band offset
MR      = iparams{1}(5);        % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
INT     = iparams{1}(6:6+nSTATES-1);             % scalar of the peak intensity of PE curve.
FWHM    = iparams{1}(6+nSTATES:6+2*nSTATES-1);	 % scalar of the FWHM of the PE curve.
% -- Extracting the limits
lbINT   = iparams{2}(6:6+nSTATES-1);
ubINT   = iparams{3}(6:6+nSTATES-1);
lbBOFF   = iparams{2}(4);
ubBOFF   = iparams{3}(4);
% -- Validity check on the inputs
if MR < 0; MR = 0; end
if MR > 1; MR = 1; end
if FDT < 0; FDT = 0; end
if FDW < 0; FDW = 0; end
for i = 1:length(INT)
    if INT(i) < 0; INT(i) = 0; end
    if lbINT(i) < 0; lbINT(i) = 0; end
    if ubINT(i) < 0; ubINT(i) = 0; end
end

%% - 1 - Extracting all background subtraction data and information
% -- Extracting the background
[X, D, B]  = PESBackground(arpesStr.xdat, arpesStr.ydat,...
    bTYPE, ibgrnd{1}(1), ibgrnd{1}(2), ibgrnd{1}(3), ibgrnd{4});

%% - 2 - Extracting all model data and information
% -- Extracting the original XPS data
xdat   	= arpesStr.xdat;
ydat 	= arpesStr.ydat;
% -- Initialising the vectors to contain XPS data
DB    = D - B;                  % Data after being background subtracted
M     = zeros(size(X));      	% XPS model data
% -- Filing through all curve components and extracting them seperately
cYY = []; BE = [];
for i = 1:nSTATES
    BE(i) = modelFunc{i}(BOFF);
    cYY(:,i) = PESCurve_FDD(X, cTYPE(i), BE(i), INT(i), FWHM(i), MR, 0, 0, 0, 0, FDEF, FDT, FDW);
    M = M + cYY(:,i);
    lbBE(i) = modelFunc{i}(lbBOFF);
    ubBE(i) = modelFunc{i}(ubBOFF);
end

%% - 3 - Determination of the residuals and chi-squared
R           = M - (D - B);	% Residuals
CHISQ       = sum(R.^2 ./ abs(M));          % Chi-squared

%% - 4 - Plotting the model to be used for fitting
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'Initial PES Curve Fitting');
fig.Position(3) = 2.5*pp.fig5x4(1); 
fig.Position(4) = pp.fig5x4(2);
%% - 4.1 - PLOTTING THE RAW DATA AND BEST FIT BACKGROUND
subplot(121); hold on;
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
subplot(4,2,[2,4,6]); hold on;
% -- Plotting all of the curve components
for i = 1:nSTATES
    area(X, cYY(:,i), 'FaceColor', pp.col.fit{i}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
    plot(X, cYY(:,i), 'k-', 'linewidth', 0.25);
end
% -- Plotting the curve component energy locations
for i = 1:nSTATES
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    line([1,1]*BE(i), [0, max(cYY(:,i))],...
        'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '-');
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
title("Initial Curves for Fitting", 'interpreter', 'none', 'fontsize', 9);
%% - 4.3 - PLOTTING THE RESIDUALS FOR THE BEST FIT TO THE DATA
subplot(4,2,8); hold on;
bar(X, R, 'facecolor', [0 0 0]);
gca_props(); grid on;
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');

%% - 5 - Plotting the uncertainties to allow for easier optimisation
subplot(4,2,[2,4,6]); hold on;
for i = 1:nSTATES
    %% - 5.1 - Plotting the primary peak uncertainties
    x_vals = [lbBE(i), lbBE(i), ubBE(i), ubBE(i), lbBE(i)];
    y_vals = [lbINT(i), ubINT(i), ubINT(i), lbINT(i), lbINT(i)];
    patch(x_vals, y_vals, pp.col.fit{i}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);
    %% - 5.2 - Adding text to identify each peak index
    text(max(x_vals), max(y_vals), string(i), 'interpreter', 'latex', 'fontsize', 14, 'color', 'k');
end

end