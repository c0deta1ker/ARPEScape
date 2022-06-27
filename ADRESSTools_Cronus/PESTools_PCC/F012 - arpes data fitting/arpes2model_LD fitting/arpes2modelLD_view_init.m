function fig = arpes2model_view_init(arpesStr, modelStr, iparams, ibgrnd, kx_lims, eb_lims)
% fig = arpes2model_view_init(arpesStr, modelStr, iparams, ibgrnd, kx_lims, eb_lims)
%   This function is used to plot the initial curve fitting model prior to
%   using the 'arpes2model_solver()' algorithm. This is used as an informative 
%   plot that allows you to view and create a better initial guess of the
%   model prior to running the fitting algorithm.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   arpesStr:      	MATLAB data-structure that contains the ARPES data.
%   -   modelStr:      	MATLAB data-structure that contains the model data.
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x6 array: the model fit parameters [K0,KFWHM,EFWHM,MFP,BOFF,INT]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x3 vector of the background parameters: [MX,MY,C]
%   -   kx_lims:      	[1x2] vector of [minKx, maxKx], which defines the consistent kx fit window
%   -   eb_lims:        [1x2] vector of [minEb, maxEb], which defines the consistent eb fit window
%
%   OUT:
%   -   fig:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters
if nargin < 5; kx_lims = 0.95*[min(modelStr.kx(:)), max(modelStr.kx(:))]; end
if nargin < 6; eb_lims = 0.95*[min(modelStr.eb(:)), max(modelStr.eb(:))]; end
if isempty(kx_lims); kx_lims = 0.95*[min(modelStr.kx(:)), max(modelStr.kx(:))]; end
if isempty(eb_lims); eb_lims = 0.95*[min(modelStr.eb(:)), max(modelStr.eb(:))]; end
% -- Validity check on the input variables
if kx_lims(1) < 0.95*min(modelStr.kx(:)); kx_lims(1) = 0.95*min(modelStr.kx(:)); end
if kx_lims(2) > 0.95*max(modelStr.kx(:)); kx_lims(2) = 0.95*max(modelStr.kx(:)); end
if eb_lims(1) < 0.95*min(modelStr.eb(:)); eb_lims(1) = 0.95*min(modelStr.eb(:)); end
if eb_lims(2) > 0.95*max(modelStr.eb(:)); eb_lims(2) = 0.95*max(modelStr.eb(:)); end

%% Initialising variables
K0      = iparams{1}(1);     % scalar of kx-shift of the ARPES data
MFP     = iparams{1}(4);     % scalar of mean free path (associated with resolution of bands)
BOFF    = iparams{1}(5);     % scalar of band-offset
% -- Defining the initial conditions and limits of the parameters
x0      = [iparams{1}(:); ibgrnd{1}(:)];
lb      = [iparams{2}(:); ibgrnd{2}(:)];
ub      = [iparams{3}(:); ibgrnd{3}(:)];

%% 1 - Making the ARPES and MODEL data consistent with the fitting window
[ARPES, MODEL] = extract_arpes_and_model(x0, arpesStr, modelStr, kx_lims, eb_lims);

%% - 1 - Determination of the residuals and chi-squared (the minimisation variable)
% - 1.1 - Extracting the DATA
D       = ARPES.data;
D       = DifC(D, 1); D = D ./ max(D(:));
% - 1.2 - Extracting the MODEL
M       = model(x0, MODEL);
M       = DifC(M, 1); M = M ./ max(M(:));
% - 1.3 - Extracting the BACKGROUND
B       = background(x0, MODEL);
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
win     = 0.01;
val1    = 0;
val2    = 0.25*min(ARPES.kx(:));
val3    = 0.25*max(ARPES.kx(:));
pp      = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'ARPES2MODEL Fitting');
fig.Position(3) = 4.*pp.fig5x4(1);
fig.Position(4) = 2.*pp.fig5x4(2);
%% - 2.1 - Plotting the model data
% -- Extracting the data matrix
data            = M + B;
% -- Plotting the data
axs(1) = subplot(3,7,[1,2,8,9]); hold on;
ImData(MODEL.kx, MODEL.eb, data); 
img_props(); 
xlabel(' $$ k_{//} [A^{-1}] $$', 'interpreter', 'latex');
ylabel(' $$ E - E_F [eV] $$', 'interpreter', 'latex');
axis([MODEL.kx_lims, MODEL.eb_lims]);
caxis([min(data(:)), max(data(:))]);
text(0.04, 0.94, "$$ \Phi = $$" + string(round(MODEL.BOFF,3)) + " $$ eV $$",...
    'interpreter', 'latex', 'fontsize', 16, 'color', 'w', 'Units','normalized');
axis square;
title('Initial model data');

%% - 2.2 - Plotting the ARPES data
axs(2) = subplot(3,7,[3,4,10,11]); hold on;
ImData(ARPES.kx, ARPES.eb, D); img_props(); 
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
[ARPES_XCut, ARPES_DCut] = Cut(ARPES.kx, ARPES.eb, D, 'edc', val1 + win*[-1, 1]);
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
[ARPES_XCut, ARPES_DCut] = Cut(ARPES.kx, ARPES.eb, D, 'edc', val2 + win*[-1, 1]);
[MODEL_XCut, MODEL_DCut] = Cut(MODEL.kx, MODEL.eb, (M + B), 'edc', val2 + win*[-1, 1]);
% -- Plotting the experimental and fit spectra
plot(ARPES_XCut, ARPES_DCut, 'k-', 'color', pp.col.fit{2}, 'linewidth', 2*pp.llwidth);
plot(MODEL_XCut, MODEL_DCut, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Extracting the RHS EDCs
[ARPES_XCut, ARPES_DCut] = Cut(ARPES.kx, ARPES.eb, D, 'edc', val3 + win*[-1, 1]);
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
    kxFWHM  = x(2);    % scalar of kx-resolution (Gaussian FWHM of k-broadening)
    ebFWHM  = x(3);    % scalar of eb-resolution (Gaussian FWHM of eb-broadening)
    MFP     = x(4);    % scalar of mean free path (associated with resolution of bands)
    BOFF    = x(5);    % scalar of band-offset
    INT     = x(6);    % scalar of the intensity scaling factor
    % 2 - Extracting the model curve based on the MFP and BOFF chosen
    model_data   	= MODEL.data;
    % 3 - Gaussian broadening the model in x- and y-dimensions
    % -- Extracting the FWHM of the smoothing
    kxFWHM         = kxFWHM ./ abs(MODEL.kx_step);
    ebFWHM         = ebFWHM ./ abs(MODEL.eb_step);
    % -- Gaussian smoothing
    model_data   	= Gaco2(model_data, kxFWHM, ebFWHM);
    % 4 - Intensity scaling
    model_data      = INT .* model_data;
    % Assigning the final model
    M               = model_data;
    M(isnan(M))     = 0;
end

%% DEFINING THE FUNCTION THAT DETERMINES THE PLANAR BACKGROUND
function B = background(x, MODEL)
    % 1 - Initialising variables
    MX  = x(end-2);   	% scalar of the x-gradient of the plane.
    MY  = x(end-1);    	% scalar of the y-gradient of the plane.
    C   = x(end);     	% scalar of the constant off-set of the plane.
    % 2 - Creating a planar background to be subtracted
    B  	= MX * MODEL.kx + MY * MODEL.eb + C;
    B(isnan(B)) = 0;
end

%% DEFINING THE FUNCTION THAT EXTRACTS THE ARPES AND MODEL DATA CONSISTENTLY
function [ARPES, MODEL] = extract_arpes_and_model(x, arpesStr, modelStr, kx_lims, eb_lims)
    % 1 - Initialising variables
    K0      = x(1);    % scalar of kx-resolution (Gaussian FWHM of k-broadening)
    MFP     = x(4);    % scalar of mean free path (associated with resolution of bands)
    BOFF    = x(5);    % scalar of band-offset
    % 2 - Cropping the ARPES data to extract dimensional size
    [xDat_crop, yDat_crop, ~] = data_crop2D(arpesStr.kx, arpesStr.eb, arpesStr.data, kx_lims, eb_lims);
    kx_length   = size(xDat_crop,2);
    eb_length   = size(yDat_crop,1);
    % 3 - Extracting the interpolated MODEL data
    imodelStr	= mstheory_interp(modelStr, BOFF, MFP);
    % 4 - Making the ARPES and MODEL data have a consistent domain
    KX          = linspace(kx_lims(1), kx_lims(2), kx_length);
    EB          = linspace(eb_lims(1), eb_lims(2), eb_length)';
    ARPES_DATA  = interpn(arpesStr.eb(:,1), arpesStr.kx(1,:)-K0, arpesStr.data, EB, KX);
    MODEL_DATA  = interpn(imodelStr.eb, imodelStr.kx, imodelStr.data, EB, KX);
    % 5 - Data renormalisation
    ARPES_DATA = ARPES_DATA - min(ARPES_DATA(:)); ARPES_DATA = ARPES_DATA ./ max(ARPES_DATA(:));
    MODEL_DATA = MODEL_DATA - min(MODEL_DATA(:)); MODEL_DATA = MODEL_DATA ./ max(MODEL_DATA(:));
    ARPES_DATA(isnan(ARPES_DATA)) = 0;
    MODEL_DATA(isnan(MODEL_DATA)) = 0;
    % 6 - Defining the ARPES data structure
    ARPES.kx           	= KX;
    ARPES.kx_step       = mean(diff(KX(:)));
    ARPES.kx_lims      	= kx_lims;
    ARPES.eb           	= EB;
    ARPES.eb_step       = mean(diff(EB(:)));
    ARPES.eb_lims     	= eb_lims;
    ARPES.data         	= ARPES_DATA;
    % 7 - Defining the MODEL data structure
    MODEL.MFP          	= MFP;
    MODEL.BOFF        	= BOFF;
    MODEL.kx           	= KX;
    MODEL.kx_step      	= mean(diff(KX(:)));
    MODEL.kx_lims      	= kx_lims;
    MODEL.eb           	= EB;
    MODEL.eb_step     	= mean(diff(EB(:)));
    MODEL.eb_lims     	= eb_lims;
    MODEL.data        	= MODEL_DATA;
end
