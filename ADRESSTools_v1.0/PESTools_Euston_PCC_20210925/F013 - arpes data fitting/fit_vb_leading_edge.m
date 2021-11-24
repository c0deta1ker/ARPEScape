function fitStr = fit_vb_leading_edge(XDat, YDat, DDat, vb_fit_args, plot_results)
% fitStr = xps_solver(xpsStr, cTYPE, iPESCurves, dPESCurves, bTYPE, iPESBgrnd, dPESBgrnd, solve_type)
%   Function that runs a a Global Optimisation Algorithm based on Simulated
%   Annealing to determine the best curve fit to the input XPS data and
%   model (initial conditions and bounds). The temperature of the Simulated
%   Annealing is set to walk a large span of the parameter space to ensure
%   the global minimum is obtained. The minimisation function is based on
%   minimising the standard deviation of the sum of the squared residuals.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   XDat:           MATLAB data-structure that contains the XPS data.
%   -   YDat:           1xN vector of the type of curve to use for the nth state.
%   -   DDat:           Nx8 array: the n'th peak parameters [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   vb_fit_args: 	Nx8 array: the n'th peak parameters uncertainties [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   plot_results: 	string of the type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters
if nargin < 5; plot_results = 1; end
if nargin < 4; vb_fit_args = {[-0.05, 0.05], [-0.25, 0.25], -0.75 + [-0.25, 0.25], []}; end
if isempty(plot_results); plot_results = 1; end
if isempty(vb_fit_args); vb_fit_args = {[-0.05, 0.05], [-0.25, 0.25], -0.75 + [-0.25, 0.25], []}; end

%% - 1 - Initialising the fitting parameters
pp          = plot_props();
% - Extracting the VB fit parameters
kWin_edc    = vb_fit_args{1}; if isempty(kWin_edc); kWin_edc = [-0.05, 0.05]; end               % window of the EDC taken around the VBM position
eWin_back 	= vb_fit_args{2}; if isempty(eWin_back); eWin_back = [-0.25, 0.25]; end             % window of the background region within the EDC cut
eWin_edge 	= vb_fit_args{3}; if isempty(eWin_edge); eWin_edge = -0.75 + [-0.25, 0.25]; end     % window of the VBM leading edge region within the EDC cut
dESmooth 	= vb_fit_args{4}; if isempty(dESmooth); dESmooth = []; end                          % can be empty []. Single, constant value that defines Gaussian pre-smoothing energy width.
% - Validity check on the inputs
kWin_edc    = sort(kWin_edc);
eWin_back   = sort(eWin_back);
eWin_edge 	= sort(eWin_edge);

%% - 2 - Extracting the EDC cut from the data
% - Extracting the full EDC cut through the VBM
[XCut, DCut] = Cut(XDat, YDat, DDat, 'edc', kWin_edc);
DCut = DCut ./ max(DCut(:));
% - Extracting the region of the background
[~, lb_indx] = min(abs(XCut - eWin_back(1)));
[~, ub_indx] = min(abs(XCut - eWin_back(2)));
XCut_back = XCut(lb_indx:ub_indx);
DCut_back = DCut(lb_indx:ub_indx);
% - Extracting the region of the edge
[~, lb_indx] = min(abs(XCut - eWin_edge(1)));
[~, ub_indx] = min(abs(XCut - eWin_edge(2)));
XCut_edge = XCut(lb_indx:ub_indx);
DCut_edge = DCut(lb_indx:ub_indx);

%% - 3 - Fitting a linear equation to the background and edge
XX          = linspace(min(XCut(:)), max(XCut(:)), 1e3)';
% -- Fitting the data
fit_back    = fit(XCut_back, DCut_back, 'poly1');
fit_edge    = fit(XCut_edge, DCut_edge, 'poly1');
% - Extracting confidence interval for each coefficient
fitci_back 	= abs(0.5*range(confint(fit_back)));
fitci_edge 	= abs(0.5*range(confint(fit_edge)));
% - Evaluating the best fit lines over a consistent domain
DD_back     = fit_back(XX);
DD_edge     = fit_edge(XX);
% - Extracting confidence interval to be plotted as a patch
ci_back     = predint(fit_back, XX, 0.95, 'observation', 'off');
ci_edge     = predint(fit_edge, XX, 0.95, 'observation', 'off');
ciXX        = [XX;         flipud(XX)];
ciDD_back	= [ci_back(:,1); flipud(ci_back(:,2))];
ciDD_edge	= [ci_edge(:,1); flipud(ci_edge(:,2))];

%% - 4 - Finding the point of intersection and its uncertainty
% - Finding the mean value of the POI
a       = fit_back.p1;
c       = fit_back.p2;
b       = fit_edge.p1;
d       = fit_edge.p2;
coeff   = (d-c)./(a-b);
X0      = round(coeff, 3);
Y0      = round(a.*coeff + c, 3);
% - Finding the uncertainty in the POI
A       = fit_back.p1 + fitci_back(1);
C       = fit_back.p2 - fitci_back(2);
B       = fit_edge.p1 - fitci_edge(1);
D       = fit_edge.p2 + fitci_edge(2);
coeff   = (D-C)./(A-B);
dX0     = round(abs(X0 - (coeff)), 3);
dY0     = round(abs(Y0 - (A.*coeff + C)), 3);

%% - 5 - Plotting the fitting result
if plot_results == 1
    % -- Initialising the figure
    fig = figure(); 
    fig.Position(3) = 1.25*pp.fig5x4(1); 
    fig.Position(4) = 2.00*pp.fig5x4(2);
    
    % -- Plotting the ARPES data and the EDC cut location
    subplot(2,1,1); hold on;
    ImData(YDat, XDat, DDat);
    patch([-1e3, 1e3, 1e3, -1e3, -1e3],[kWin_edc(1), kWin_edc(1), kWin_edc(2), kWin_edc(2), kWin_edc(1)],...
        [0 1 0], 'linewidth', 1, 'facealpha', 0, 'edgecolor', [0 0 1]);
    img_props();
    axis([min(YDat(:)), max(YDat(:)),min(XDat(:)), max(XDat(:))]);
    xlabel('$$ \bf  E_{B} (eV) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
    
    % -- Plotting the EDC cut and the VB fitting to the leading edge 
    subplot(2,1,2); hold on;
    % -- Plotting the full EDC cut first
    plot(XCut, DCut, 'k.-', 'color', pp.col.fit{1}, 'LineWidth', 0.50);
    % -- Plotting the fits to the EDC background
    patch(ciXX, ciDD_back, pp.col.fit{2}, 'EdgeAlpha', 0, 'facealpha', 0.25);
    plot(XCut_back, DCut_back, 'k.-', 'color', pp.col.fit{2}, 'LineWidth', 1.50);
    plot(XX, DD_back, 'k:', 'color', pp.col.fit{2}, 'LineWidth', 1.50);
    % -- Plotting the fits to the EDC leading edge
    patch(ciXX, ciDD_edge, pp.col.fit{3}, 'EdgeAlpha', 0, 'facealpha', 0.25);
    plot(XCut_edge, DCut_edge, 'g.-', 'color', pp.col.fit{3}, 'LineWidth', 1.50);
    plot(XX, DD_edge, 'g:', 'color', pp.col.fit{3}, 'LineWidth', 1.50);
    % -- Plotting the point of intersection    
    errorbar(X0, Y0, dY0, dY0, dX0, dX0, 'ko',...
        'markersize', 7, 'color', [0 0 0], 'markerfacecolor', [1 0 0]);
    % -- Add text to show the POI
    text(0.05, 0.90, "$$ E_{VBM} = (" + X0 +  " \pm " + dX0 + ") eV $$",...
        'interpreter', 'latex', 'fontsize', 14, 'color', 'k', 'Units','normalized');
    % -- Formatting the figure
    gca_props();
    ax = gca;
    ax.XAxisLocation = 'bottom';            % 'bottom' | 'top' | 'origin'
    ax.YAxisLocation = 'right';             % 'left' | 'right' | 'origin'
    % - Axis labels and limits
    xlabel('$$ \bf  E_{B} (eV) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
    % - Plotting the x- and y-axes
    axis([min(YDat(:)), max(YDat(:)), 0, 1.15]);
end


%% 6 - Appending data to MATLAB data structure
fitStr              = struct();
fitStr.vb_fit_args  = vb_fit_args;
% - Appending all the EDC cuts taken
fitStr.XCut         = XCut;
fitStr.DCut         = DCut;
fitStr.XCut_back  	= XCut_back;
fitStr.DCut_back  	= DCut_back;
fitStr.XCut_edge  	= XCut_edge;
fitStr.DCut_edge  	= DCut_edge;
% - Appending the best fits to the background and edge
fitStr.fit_back  	= fit_back;
fitStr.fit_edge  	= fit_edge;
fitStr.XX           = XX;
fitStr.DD_back  	= DD_back;
fitStr.DD_edge  	= DD_edge;
% - Appending the VBM position and its uncertainty
fitStr.VBM          = X0;
fitStr.dVBM         = dX0;
fitStr.VBM_int      = Y0;
fitStr.dVBM_int     = dY0;

end