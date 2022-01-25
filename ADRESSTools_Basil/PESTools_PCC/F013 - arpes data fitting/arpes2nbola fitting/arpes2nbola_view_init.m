function fig = arpes2nbola_view_init(arpesStr, cTYPE, iparams, ibgrnd)
% fig = arpes2nbola_view_init(arpesStr, cTYPE, iparams, ibgrnd)
%   This function is used to plot the initial guess of the ARPES curve fitting
%   performed by 'arpes2nbola_solver()'. The plot consists of 3 subplots; (1) The
%   initial model data to be fitted based on the initial parmaeters; (2) A 
%   plot showing all of the experimental ARPES data to be fitted; (3) A 
%   plot summarising the residuals, showing the quality of the experimental
%   and model fit.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   arpesStr:       MATLAB data-structure that contains the initial ARPES data.
%   -   cTYPE:          1xN vector of the type of curve to use for fitting. Default: "G2DA" ("G2D", "L2D")
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x6 array: the model fit parameters [INT,XLOC,YLOC,XFWHM,YFWHM,MSTAR]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x6 vector of the background parameters: [FDEF,FDT,FDW,BGR,BIN,BCO]
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Initialising variables
% - Defining the initial conditions and limits of the parameters
x0      = [iparams{1}(:); ibgrnd{1}(:)];
lb      = [iparams{2}(:); ibgrnd{2}(:)];
ub      = [iparams{3}(:); ibgrnd{3}(:)];
% - Making the ARPES and MODEL data consistent over a consistent domain
[ARPES, MODEL] = extract_arpes_and_model(x0, arpesStr, cTYPE);
% - Extracting the parabolic dispersions for each subband energy
sbn_kx      = linspace(ARPES.kx_lims(1), ARPES.kx_lims(2), 1e3);
sbn_eb      = {};
XLOC    = iparams{1}(:,2);
YLOC    = iparams{1}(:,3);
MSTAR   = iparams{1}(:,6);
for i = 1:length(YLOC)
    sbn_eb{i}      = eb_calc(sbn_kx, XLOC(i), YLOC(i), MSTAR(i));
end

%% - 1 - Determination of the residuals and chi-squared (the minimisation variable)
% - 1.1 - Extracting the DATA
D       = ARPES.data;
% - 1.2 - Extracting the MODEL
M       = MODEL.data;
% - 1.4 - Extracting the RESIDUALS
R       = D - M;            % residuals = data - model
Rx      = mean(R, 1);       % R in 1D along x
Ry  	= mean(R, 2);       % R in 1D along y
% - 1.5 - Extracting the CHI-SQUARED value to be minimised
CHISQ2D     = R.^2 ./ M;
CHISQxx   	= mean(CHISQ2D, 1);      % CHISQ in 1D along x
CHISQyy  	= mean(CHISQ2D, 2);      % CHISQ in 1D along y
% - 1.6 - Extracting the REDUCED CHI-SQUARED value including degrees of freedom
CHISQ       = sum(sum(CHISQ2D));

%% - 2 - PLOTTING THE FINAL DATA COMPARED WITH MODEL AT BOFF = 0 eV
% -- Initialising the plot properties
win     = 0.01;
val1    = 0;
val2    = 0.15*min(ARPES.kx(:));
val3    = 0.15*max(ARPES.kx(:));
pp      = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'Initial Model');
fig.Position(3) = 4.*pp.fig5x4(1);
fig.Position(4) = 2.*pp.fig5x4(2);
%% - 2.1 - Plotting the model data
% -- Extracting the data matrix
data            = M;
% -- Plotting the data
axs(1) = subplot(3,7,[1,2,8,9]); hold on;
ImData(MODEL.kx, MODEL.eb, data); 
for i = 1:length(sbn_eb); plot(sbn_kx, sbn_eb{i}, 'r:', 'linewidth', 3); end
img_props(); 
xlabel(' $$ k_{//} [A^{-1}] $$', 'interpreter', 'latex');
ylabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
axis([MODEL.kx_lims, MODEL.eb_lims]);
caxis([min(data(:)), max(data(:))]);
axis square;
title('Initial model data');

%% - 2.2 - Plotting the ARPES data
axs(2) = subplot(3,7,[3,4,10,11]); hold on;
ImData(ARPES.kx, ARPES.eb, ARPES.data); img_props(); 
for i = 1:length(sbn_eb); plot(sbn_kx, sbn_eb{i}, 'r:', 'linewidth', 3); end
xlabel(' $$ k_{//} [A^{-1}] $$', 'interpreter', 'latex');
ylabel('', 'interpreter', 'latex');
axis([ARPES.kx_lims, ARPES.eb_lims]); axis square;
caxis([min(ARPES.data(:)), max(ARPES.data(:))]);
title('ARPES data to be fitted');
% -- Plotting the EDC cut lines
line([1 1].*val1, [-1e5, 1e5], 'Color', pp.col.fit{1}, 'LineWidth', pp.llwidth, 'Linestyle', '-');
line([1 1].*val2, [-1e5, 1e5], 'Color', pp.col.fit{2}, 'LineWidth', pp.llwidth, 'Linestyle', '-');
line([1 1].*val3, [-1e5, 1e5], 'Color', pp.col.fit{3}, 'LineWidth', pp.llwidth, 'Linestyle', '-');

%% - 2.3 - Plotting the 2D residuals
axs(3) = subplot(3,7,[5,6,12,13]); hold on;
ImData(ARPES.kx, ARPES.eb, CHISQ2D); img_props(); 
axis([ARPES.kx_lims, ARPES.eb_lims]);
title('Chi-Squared');
% -- Add annotation for the quality of fit
text(0.04, 0.94, "$$ \chi^2 = $$ " + string(CHISQ),...
    'interpreter', 'latex', 'fontsize', 16, 'color', 'w', 'Units','normalized');
% -- Make the colormap bipolar (negative values blue, positive is red)
colormap(axs(3),'jet');
% -- Make the colormap bipolar (negative values blue, positive is red)
ylabel('', 'interpreter', 'latex');
xlabel('');

%% - 2.4 - Plotting the 1D residuals along x
axs(4) = subplot(3,7,[19,20]); hold on;
bar(ARPES.kx, CHISQxx);
gca_props(); grid on;
xlabel(' $$ k_{//} [A^{-1}] $$', 'interpreter', 'latex');
ylabel('$$ \bf  \chi^2(x) $$', 'Interpreter', 'latex');
xlim(ARPES.kx_lims);
ax = gca; ax.YAxisLocation = 'right';

%% - 2.5 - Plotting the 1D residuals along y
axs(5) = subplot(3,7,[7,14]); hold on;
barh(ARPES.eb, CHISQyy);
gca_props(); grid on;
ylabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
xlabel('$$ \bf  \chi^2(y) $$', 'Interpreter', 'latex');
ylim(ARPES.eb_lims);
ax = gca; ax.YAxisLocation = 'right';

%% - 2.6 - Plotting best fit EDCs
%% - 2.6.1 - Gamma EDCs
axs(6) = subplot(3,7,[15,16]); hold on;
% -- Extracting the EDCs
[ARPES_XCut, ARPES_DCut] = Cut(ARPES.kx, ARPES.eb, ARPES.data, 'edc', val1 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(MODEL.kx, MODEL.eb, MODEL.data, 'edc', val1 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{1}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Plotting the sub-band energies
for i = 1:length(YLOC); line([1 1]*YLOC(i), [-1e5, 1e5], 'Color', 'r', 'LineWidth', 2, 'Linestyle', ':'); end
gca_props(); grid on;
ylabel('$$ Intensity $$', 'Interpreter', 'latex');
xlabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
axis([min(ARPES_XCut(:)), max(ARPES_XCut(:)), 0, 1.25*max(ARPES_DCut(:))]);
%% - 2.6.1 - LHS and RHS EDCs
axs(7) = subplot(3,7,[17,18]); hold on;
% -- Extracting the LHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(ARPES.kx, ARPES.eb, ARPES.data, 'edc', val2 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(MODEL.kx, MODEL.eb, MODEL.data, 'edc', val2 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{2}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Extracting the RHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(ARPES.kx, ARPES.eb, ARPES.data, 'edc', val3 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(MODEL.kx, MODEL.eb, MODEL.data, 'edc', val3 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{3}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
gca_props(); grid on;
xlabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
axis([min(ARPES_XCut(:)), max(ARPES_XCut(:)), 0, 1.25*max(ARPES_DCut(:))]);

end

%% DEFINING THE FUNCTION THAT EXTRACTS THE ARPES AND MODEL DATA CONSISTENTLY
function [ARPES, MODEL] = extract_arpes_and_model(x, arpesStr, cTYPE)
    % 1 - Initialising variables
    n = length(cTYPE);
    % - Extracting parabolic parameters
    INT     = x(1:n);           % scalar of the peak intensity of 2D PE curve.
    XLOC    = x(n+1:2*n);       % scalar of the x-location of the min/max of the 2D parabolic ARPES dispersion.
    YLOC    = x(2*n+1:3*n);  	% scalar of the y-location of the min/max of the 2D parabolic ARPES dispersion.
    XFWHM   = x(3*n+1:4*n);     % scalar of the x-axis FWHM for each Gaussian (k-resolution)
    YFWHM   = x(4*n+1:5*n);   	% scalar of the y-axis FWHM for each Gaussian (Eb-resolution)
    MSTAR   = x(5*n+1:6*n);   	% scalar of the effective mass, which governs the curvature of the parabola.
    % - Extracting background parameters
    FDEF    = x(end-5);         % scalar of the FDD Fermi-Level position.
    FDT     = x(end-4);         % scalar of the FDD temperature.
    FDW     = x(end-3);         % scalar of the FDD Gaussian width after convolution.
    BGR     = x(end-2);         % scalar of the gradient of the linear background.
    BIN     = x(end-1);         % scalar of the y-intercept of the linear background.
    BCO     = x(end);           % scalar of the constant background y-offset value.
    % - Extracting the kx and eb limits
    kx_lims = [min(arpesStr.kx(:)), max(arpesStr.kx(:))];
    eb_lims = [min(arpesStr.eb(:)), max(arpesStr.eb(:))];
    % 2 - Extracting the constructed MODEL data
    KX    = linspace(kx_lims(1), kx_lims(2), size(arpesStr.kx,2));
    EB    = linspace(eb_lims(1), eb_lims(2), size(arpesStr.eb,1))';
    MODEL_DATA    = ARPESCurve_FDDGpL(KX, EB, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR, FDEF, FDT, FDW, BGR, BIN, BCO);
    ARPES_DATA    = arpesStr.data;
    % - Removing any NaN values
    MODEL_DATA(isnan(MODEL_DATA)) = 0;
    ARPES_DATA(isnan(ARPES_DATA)) = 0;
    % 3 - Defining the ARPES data structure
    ARPES               = struct();
    ARPES.kx           	= KX;
    ARPES.kx_step       = mean(diff(KX(:)));
    ARPES.kx_lims      	= kx_lims;
    ARPES.eb           	= EB;
    ARPES.eb_step       = mean(diff(EB(:)));
    ARPES.eb_lims     	= eb_lims;
    ARPES.data         	= ARPES_DATA;
    % 4 - Defining the MODEL data structure
    MODEL               = struct();
    MODEL.kx           	= KX;
    MODEL.kx_step      	= mean(diff(KX(:)));
    MODEL.kx_lims      	= kx_lims;
    MODEL.eb           	= EB;
    MODEL.eb_step     	= mean(diff(EB(:)));
    MODEL.eb_lims     	= eb_lims;
    MODEL.data        	= MODEL_DATA;
end