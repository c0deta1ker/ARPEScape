function fig = curvefit1D_view_init(xdat, ydat, curve_type, init_params, bgrnd_type, init_bgrnd)
% fig = curvefit1D_view_init(xdat, ydat, curve_type, init_params, bgrnd_type, init_bgrnd)
%   This function is used to plot the initial curve fitting model prior to
%   using the 'curvefit1D_solver()' algorithm. This is used as an informative 
%   plot that allows you to view and create a better initial guess of the
%   model prior to running the fitting algorithm.
%   
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   xdat:           [nX x 1] array of the x-axis data
%   -   ydat:           [nY x 1] array of the y-axis data
%   -   curve_type:  	string of the type of curve to be used for fitting. Default: "sGLA" ("pGLA", "sGL", "pGL", "G", "L", "DS")
%   -   init_params:   	{3 x 1} cell array of {x0}{lb}{ub}, where each cell is an [N x 8] matrix of N'th peak parameters: [XLOC,YINT,FWHM,MR,ASY,LSX,LSY,LSW]
%   -   [bgrnd_type]:  	string of the type of background to use for fitting. Default: "none" ("Poly", "Shir", "LinShir")
%   -   [init_bgrnd]:   {3 x 1} cell array of {x0}{lb}{ub}, where each cell is an [1 x 5] matrix of the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%
%   OUT:
%   -   fig:            MATLAB figure object showing the initial model

%% Default input parameters & consistency checks
pp  = plot_props();
% -- Defining input parameters if undefined
if nargin < 6; init_bgrnd{1} = [min(xdat(:)), max(xdat(:)), 1, 0.5, 0, 0]; init_bgrnd{2} = init_bgrnd{1}; init_bgrnd{2}(6) = min(ydat(:)); init_bgrnd{3} = init_bgrnd{1}; init_bgrnd{3}(6) = max(ydat(:)); end
if nargin < 5; bgrnd_type = "none"; end
if isempty(init_bgrnd); init_bgrnd{1} = [min(xdat(:)), max(xdat(:)), 1, 0.5, 0, 0]; init_bgrnd{2} = init_bgrnd{1}; init_bgrnd{2}(6) = min(ydat(:)); init_bgrnd{3} = init_bgrnd{1}; init_bgrnd{3}(6) = max(ydat(:)); end
if isempty(bgrnd_type); bgrnd_type = "none"; end
if isempty(curve_type); curve_type = "sGLA"; end
% -- Consistency check and finding the total number of curves
if size(init_params{1}, 1) ~= size(init_params{2}, 1) || size(init_params{2}, 1) ~= size(init_params{3}, 1)
    error('The input parameter cell array is not a consistent size - check iparams input!');
end
% -- Extracting the total number of curves
ncurves = size(init_params{1}, 1);
% -- Validity check on input parameters
for i = 1:ncurves
    if size(init_params{1}, 2) == 4
        % --- Set the ASY, LSE, LSI and LSW components to zero
        init_params{1}(i,5) = 0; init_params{1}(i,6) = 0; init_params{1}(i,7) = 0; init_params{1}(i,8) = 0;
        init_params{2}(i,5) = 0; init_params{2}(i,6) = 0; init_params{2}(i,7) = 0; init_params{2}(i,8) = 0;
        init_params{3}(i,5) = 0; init_params{3}(i,6) = 0; init_params{3}(i,7) = 0; init_params{3}(i,8) = 0;
    elseif size(init_params{1}, 2) == 5
        % --- Set the LSE, LSI and LSW components to zero
        init_params{1}(i,6) = 0; init_params{1}(i,7) = 0; init_params{1}(i,8) = 0;
        init_params{2}(i,6) = 0; init_params{2}(i,7) = 0; init_params{2}(i,8) = 0;
        init_params{3}(i,6) = 0; init_params{3}(i,7) = 0; init_params{3}(i,8) = 0;
    elseif size(init_params{1}, 2) == 6
        % --- Set the LSI and LSW components to zero
        init_params{1}(i,7) = 0; init_params{1}(i,8) = 0;
        init_params{2}(i,7) = 0; init_params{2}(i,8) = 0;
        init_params{3}(i,7) = 0; init_params{3}(i,8) = 0;
    elseif size(init_params{1}, 2) == 7
        % --- Set the LSW components to zero
        init_params{1}(i,8) = 0;
        init_params{2}(i,8) = 0;
        init_params{3}(i,8) = 0;
    elseif size(init_params{1}, 2) > 8 || size(init_params{1}, 2) < 4 
        error('Not enough input arguments defined - check iparams input!');
    end
end

%% - 1 - Extracting all of the curve parameters
for i = 1:ncurves
    XLOC(i) = init_params{1}(i,1); lbXLOC(i) = init_params{2}(i,1); ubXLOC(i) = init_params{3}(i,1);
    YINT(i) = init_params{1}(i,2); lbYINT(i) = init_params{2}(i,2); ubYINT(i) = init_params{3}(i,2);
    FWHM(i) = init_params{1}(i,3); lbFWHM(i) = init_params{2}(i,3); ubFWHM(i) = init_params{3}(i,3);
    MR(i)   = init_params{1}(i,4); lbMR(i)   = init_params{2}(i,4); ubMR(i)   = init_params{3}(i,4);
    ASY(i)  = init_params{1}(i,5); lbASY(i)  = init_params{2}(i,5); ubASY(i)  = init_params{3}(i,5);
    LSX(i)  = init_params{1}(i,6); lbLSX(i)  = init_params{2}(i,6); ubLSX(i)  = init_params{3}(i,6);
    LSY(i)  = init_params{1}(i,7); lbLSY(i)  = init_params{2}(i,7); ubLSY(i)  = init_params{3}(i,7);
    LSW(i)  = init_params{1}(i,8); lbLSW(i)  = init_params{2}(i,8); ubLSW(i)  = init_params{3}(i,8);
end

%% - 2 - Extracting all background subtraction data and information
% -- Extracting the background
[X, D, B] = PESBackground(xdat, ydat,...
    bgrnd_type, init_bgrnd{1}(1), init_bgrnd{1}(2), init_bgrnd{1}(3), init_bgrnd{1}(4), init_bgrnd{1}(5), init_bgrnd{1}(6));

%% - 3 - Extracting all model data and information
% -- Initialising the vectors to contain XPS data
DB      = D - B;                % Data after being background subtracted
M   	= zeros(size(X));      	% XPS model data
% -- Filing through all curve components and extracting them seperately
cYY = [];
for i = 1:ncurves
   	cYY(:,i) = PESCurve(X, curve_type, XLOC(i), YINT(i), FWHM(i), MR(i), LSX(i), LSY(i), LSW(i), ASY(i));
    M = M + cYY(:,i);
end

%% - 4 - Determination of the residuals and chi-squared
R       = D - (M + B);              % Residuals
CHISQ	= sum(R.^2 ./ abs(M));     	% Chi-squared

%% - 5.1 - Plotting the initial model figure
% - Finding the optimal y-limits
minY    = min([min(cYY(:)), min(M(:)), min(D(:)), min(B(:)), min(DB(:))]);
maxY    = 1.25*max([max(cYY(:)), max(M(:)), max(D(:)), max(B(:)), max(DB(:))]);
ylims   = [minY, maxY];
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'Curve Fit 1D - Initial Model');
fig.Position(3) = 2.5*pp.fig5x4(1); 
fig.Position(4) = 1.0*pp.fig5x4(2);
% -- Plotting the raw data and background
subplot(121); hold on;
h = patch(...
    [init_bgrnd{1}(1), init_bgrnd{1}(1), init_bgrnd{1}(2), init_bgrnd{1}(2), init_bgrnd{1}(1)],...
    [-1, 1, 1, -1, -1].*1e4, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', [0 0 0]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(xdat, ydat, 'b-', 'linewidth', 0.5*pp.lwidth);
plot(X, D, 'b-', 'linewidth', pp.llwidth);
plot(X, B, 'r-', 'linewidth', pp.llwidth);
plot(X, DB, 'k-', 'linewidth', pp.llwidth);
gca_props(); grid on;
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
axis([min(xdat(:)), max(xdat(:)), ylims]);
legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
title("Background Subtraction", 'interpreter', 'none', 'fontsize', 9);
% -- Plotting the best fit curve components and the initial model
subplot(4,2,[2,4,6]); hold on;
% -- Plotting all of the curve components
for i = 1:ncurves
    area(X, cYY(:,i), 'FaceColor', pp.col.fit{i}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
    plot(X, cYY(:,i), 'k-', 'linewidth', 0.25);
end
% -- Plotting the curve component energy locations
for i = 1:ncurves
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    line([init_params{1}(i,1), init_params{1}(i,1)], [0, max(cYY(:,i))],...
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
axis([min(X(:)), max(X(:)), ylims]);
title("Initial Model", 'interpreter', 'none', 'fontsize', 9);
% -- Plotting the uncertainties to allow for easier optimisation
for i = 1:ncurves
    % --- Plotting the primary peak uncertainties
    lbBE  	= init_params{2}(i,1); ubBE      = init_params{3}(i,1);
    lbINT  	= init_params{2}(i,2); ubINT     = init_params{3}(i,2);
    x_vals = [lbBE, lbBE, ubBE, ubBE, lbBE];
    y_vals = [lbINT, ubINT, ubINT, lbINT, lbINT];
    patch(x_vals, y_vals, pp.col.fit{i}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);
    % --- Adding text to identify each peak index
    text(max(x_vals), max(y_vals), string(i), 'interpreter', 'latex', 'fontsize', 14, 'color', 'k');
    % --- Plotting the spin-orbit split (SOS) peak uncertainties
    x_vals  = XLOC(i) + [lbLSX, lbLSX, ubLSX, ubLSX, lbLSX];
    y_vals  = [lbINT*lbLSY, ubINT*ubLSY, ubINT*ubLSY, lbINT*lbLSY, lbINT*lbLSY];
    patch(x_vals, y_vals, pp.col.fit{i}, 'FaceAlpha', 0.65, 'EdgeAlpha', 0);
end

% -- Plotting the residuals of the model
subplot(4,2,8); hold on;
bar(X, R, 'facecolor', [0 0 0]);
gca_props(); grid on;
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');

end