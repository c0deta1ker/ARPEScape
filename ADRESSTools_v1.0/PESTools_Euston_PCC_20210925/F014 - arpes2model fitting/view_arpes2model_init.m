% --- General function to view the initial guess model vs ARPES data prior to fitting
function fig = view_arpes2model_init(MODEL, ARPES, iparams, ibgrnd)
% fig = view_arpes2model_init(MODEL, ARPES, iparams, ibgrnd)
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
%   -   MODEL:          MATLAB data-structure that contains the model data.
%   -   ARPES:         	MATLAB data-structure that contains the ARPES data.
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x5 array: the model fit parameters [kxFWHM,ebFWHM,MFP,BOFF,INT]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x3 vector of the background parameters: [MX,MY,C]
%
%   OUT:
%   -   fig:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Initialising input variables
MFP     = iparams{1}(3);    % scalar of mean free path (associated with resolution of bands)
BOFF    = iparams{1}(4);    % scalar of band-offset
% -- Defining the initial conditions and limits of the parameters
x0      = [iparams{1}(:); ibgrnd{1}(:)];
lb      = [iparams{2}(:); ibgrnd{2}(:)];
ub      = [iparams{3}(:); ibgrnd{3}(:)];
% -- Defining the values of the EDC cut figures
win     = 0.01;
val1    = 0; 
val2    = 0.5*min(ARPES.kx(:));
val3    = 0.5*max(ARPES.kx(:));

%% - 1 - Determination of the residuals and chi-squared (the minimisation variable)
% - 1.1 - Extracting the DATA
D   = ARPES.data;
% - 1.2 - Extracting the MODEL
M   = model(x0, MODEL);
% - 1.3 - Extracting the BACKGROUND
B   = background(x0, MODEL);
% - 1.4 - Extracting the RESIDUALS
R       = D - (M + B);  	% residuals = data - (model + background)
Rx      = mean(R, 1);       % R in 1D along x
Ry  	= mean(R, 2);       % R in 1D along y
% - 1.5 - Extracting the CHI-SQUARED value to be minimised
CHISQ2D     = R.^2 ./ (M + B);
CHISQxx   	= mean(CHISQ2D, 1);      % CHISQ in 1D along x
CHISQyy  	= mean(CHISQ2D, 2);      % CHISQ in 1D along y
% - 1.6 - Extracting the REDUCED CHI-SQUARED value including degrees of freedom
CHISQ       = sum(sum(CHISQ2D));

%% - 2 - PLOTTING THE FINAL DATA COMPARED WITH MODEL AT BOFF = 0 eV
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'ARPES2MODEL Fitting');
fig.Position(3) = 4.*pp.fig5x4(1);
fig.Position(4) = 2.*pp.fig5x4(2);
%% - 2.1 - Plotting the model data
% -- Extracting the data matrix
[~, nBOFF]      = min(abs(MODEL.BOFF - BOFF));
data            = M + B;
% -- Plotting the data
axs(1) = subplot(3,7,[1,2,8,9]); hold on;
ImData(MODEL.kx, MODEL.eb, data); 
img_props(); 
xlabel(' $$ k_{//} [A^{-1}] $$', 'interpreter', 'latex');
ylabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
axis([MODEL.kx_lims, MODEL.eb_lims]);
caxis([min(data(:)), max(data(:))]);
text(0.04, 0.94, "$$ \Phi = $$" + string(round(MODEL.BOFF(nBOFF),3)) + "$$ eV $$",...
    'interpreter', 'latex', 'fontsize', 16, 'color', 'w', 'Units','normalized');
axis square;
title('Initial model data');

%% - 2.2 - Plotting the ARPES data
axs(2) = subplot(3,7,[3,4,10,11]); hold on;
ImData(ARPES.kx, ARPES.eb, ARPES.data); img_props(); 
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
[MODEL_XCut, MODEL_DCut] = Cut(MODEL.kx, MODEL.eb, (M + B), 'edc', val1 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{1}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
gca_props(); grid on;
ylabel('$$ Intensity $$', 'Interpreter', 'latex');
xlabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
%% - 2.6.1 - LHS and RHS EDCs
axs(7) = subplot(3,7,[17,18]); hold on;
% -- Extracting the LHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(ARPES.kx, ARPES.eb, ARPES.data, 'edc', val2 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(MODEL.kx, MODEL.eb, (M + B), 'edc', val2 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{2}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Extracting the RHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(ARPES.kx, ARPES.eb, ARPES.data, 'edc', val3 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(MODEL.kx, MODEL.eb, (M + B), 'edc', val3 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{3}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
gca_props(); grid on;
xlabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');

end

%% DEFINING THE FUNCTION THAT DETERMINES THE MODEL
function M = model(x, MODEL)
    % 1 - Initialising variables
    kxFWHM  = x(1);    % scalar of kx-resolution (Gaussian FWHM of k-broadening)
    ebFWHM  = x(2);    % scalar of eb-resolution (Gaussian FWHM of eb-broadening)
    MFP     = x(3);    % scalar of mean free path (associated with resolution of bands)
    BOFF    = x(4);    % scalar of band-offset
    INT     = x(5);    % scalar of the intensity scaling factor
    % 2 - Extracting the model curve based on the MFP and BOFF chosen
    [~, nMFP]       = min(abs(MODEL.MFP - MFP));
    [~, nBOFF]      = min(abs(MODEL.BOFF - BOFF));
    model_data   	= squeeze(MODEL.data(nMFP,nBOFF,:,:))';
    % 3 - Gaussian broadening the model in x- and y-dimensions
    % -- Extracting the FWHM of the smoothing
    kxFWHM         = kxFWHM ./ abs(MODEL.kx(2) - MODEL.kx(1));
    ebFWHM         = ebFWHM ./ abs(MODEL.eb(2) - MODEL.eb(1));
    % -- Gaussian smoothing
    model_data   	= Gaco2(model_data, kxFWHM, ebFWHM);
    % 4 - Intensity scaling
    model_data      = INT .* model_data;
    % Assigning the final model
    M               = model_data;
end

%% DEFINING THE FUNCTION THAT DETERMINES THE PLANAR BACKGROUND
function B = background(x, MODEL)
    % 1 - Initialising variables
    MX  = x(end-2);   	% scalar of the x-gradient of the plane.
    MY  = x(end-1);    	% scalar of the y-gradient of the plane.
    C   = x(end);     	% scalar of the constant off-set of the plane.
    % 2 - Creating a planar background to be subtracted
    B  	= MX * MODEL.kx + MY * MODEL.eb + C;
end