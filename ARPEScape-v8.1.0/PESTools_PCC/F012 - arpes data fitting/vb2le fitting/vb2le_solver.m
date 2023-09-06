function fitStr = vb2le_solver(kx, eb, data, bgrnd_type, iparams, plot_results)
% ffitStr = vb2le_solver(kx, eb, data, bgrnd_type, iparams, plot_results)
%   Function that fits to the valence band leading edge using a 1st order
%   polynomial ('poly1') for the leading edge. The background can be
%   defined also as a 1st order polynomial ('poly1') or a constant flat
%   baseline ('flat'). The intersection between these two lines then
%   determines the VBM position based on the Leading Edge (LE) method
%   discussed in literature.
%
%   IN:
%   -   kx:             vector / array of the wave vector (x-domain) of the ARPES data.
%   -   eb:             vector / array of the binding energy (y-domain) of the ARPES data.
%   -   data:           NxM array of ARPES intensity.
%   -   bgrnd_type:   	String of either "poly0" (flat background), "poly1" (linear background), "none" (to the x-axis baseline)
%   -   iparams:      	1x6 array: the model fit parameters [KWIN,EWIN_BGRND,EWIN_EDGE,dESmooth]
%   -   plot_results: 	Either 1 or 0; plots the fits if 1.
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters
if nargin < 6; plot_results = 1; end
if isempty(plot_results); plot_results = 1; end
if isempty(bgrnd_type); bgrnd_type = "poly0"; end

%% - 1 - Initialising the fitting parameters
pp          = plot_props();
% - Extracting the VB fit parameters
KWIN        = iparams{1}; if isempty(KWIN); KWIN = [-0.05, 0.05]; end                       % window of the EDC taken around the VBM position
EWIN_BGRND 	= iparams{2}; if isempty(EWIN_BGRND); EWIN_BGRND = [-0.25, 0.25]; end         	% window of the background region within the EDC cut
EWIN_EDGE 	= iparams{3}; if isempty(EWIN_EDGE); EWIN_EDGE = -0.75 + [-0.25, 0.25]; end     % window of the VBM leading edge region within the EDC cut
dESmooth 	= iparams{4}; if isempty(dESmooth); dESmooth = 0; end                          % can be empty []. Single, constant value that defines Gaussian pre-smoothing energy width.
% - Validity check on the inputs
KWIN        = sort(KWIN);
EWIN_BGRND	= sort(EWIN_BGRND);
EWIN_EDGE 	= sort(EWIN_EDGE);
dESmooth    = abs(dESmooth);
% - Capping the k-integration window to the min/max limits
if KWIN(1) < min(kx(:)); KWIN(1) = min(kx(:));end
if KWIN(2) > max(kx(:)); KWIN(2) = max(kx(:)); end

%% - 2 - Extracting the EDC cut from the data
% - Extracting the full EDC cut through the VBM
[XCut, DCut] = Cut(kx, eb, data, 'edc', KWIN);
DCut = DCut ./ max(DCut(:));
% - Smoothing the data if necessary
if dESmooth ~= 0; DCut = Gaco1(DCut, dESmooth); end
% - Extracting the region of the background
[~, lb_indx] = min(abs(XCut - EWIN_BGRND(1)));
[~, ub_indx] = min(abs(XCut - EWIN_BGRND(2)));
XCut_back = XCut(lb_indx:ub_indx);
DCut_back = DCut(lb_indx:ub_indx);
% - Extracting the region of the edge
[~, lb_indx] = min(abs(XCut - EWIN_EDGE(1)));
[~, ub_indx] = min(abs(XCut - EWIN_EDGE(2)));
XCut_edge   = XCut(lb_indx:ub_indx);
DCut_edge   = DCut(lb_indx:ub_indx);
XX          = linspace(min(XCut(:)), max(XCut(:)), 1e3)';

%% - 3 - Fitting a linear equation to the leading edge
% -- Fitting the data
fit_edge    = fit(XCut_edge, DCut_edge, 'poly1');
% - Extracting confidence interval for each coefficient
fitci_edge 	= abs(0.5*range(confint(fit_edge)));
% - Evaluating the best fit lines over a consistent domain
DD_edge     = fit_edge(XX);
% - Extracting confidence interval to be plotted as a patch
ci_edge     = predint(fit_edge, XX, 0.95, 'observation', 'off');
ciXX        = [XX;         flipud(XX)];
ciDD_edge	= [ci_edge(:,1); flipud(ci_edge(:,2))];

%% - 4 - Fitting a to the background
if bgrnd_type == "poly0"
    % -- Fitting the data
    fit_back.p1    = 0;
    fit_back.p2    = mean(DCut_back(:));
    % - Extracting confidence interval for each coefficient
    fitci_back 	= 3*std(DCut_back(:))*[1,1];
    % - Evaluating the best fit lines over a consistent domain
    DD_back     = zeros(size(XX))+fit_back.p2;
    % - Extracting confidence interval to be plotted as a patch
    ci_back     = [DD_back-fitci_back, DD_back+fitci_back];
    ciXX        = [XX;         flipud(XX)];
    ciDD_back	= [ci_back(:,1); flipud(ci_back(:,2))];
elseif bgrnd_type == "poly1"
    % -- Fitting the data
    fit_back    = fit(XCut_back, DCut_back, 'poly1');
    % - Extracting confidence interval for each coefficient
    fitci_back 	= abs(0.5*range(confint(fit_back)));
    % - Evaluating the best fit lines over a consistent domain
    DD_back     = fit_back(XX);
    % - Extracting confidence interval to be plotted as a patch
    ci_back     = predint(fit_back, XX, 0.95, 'observation', 'off');
    ciXX        = [XX;         flipud(XX)];
    ciDD_back	= [ci_back(:,1); flipud(ci_back(:,2))];
elseif bgrnd_type == "none"
    % -- Fitting the data
    fit_back.p1    = 0;
    fit_back.p2    = 0;
    % - Extracting confidence interval for each coefficient
    fitci_back 	= 3*std(0.*DCut_back(:))*[1,1];
    % - Evaluating the best fit lines over a consistent domain
    DD_back     = zeros(size(XX))+fit_back.p2;
    % - Extracting confidence interval to be plotted as a patch
    ci_back     = [DD_back-fitci_back, DD_back+fitci_back];
    ciXX        = [XX;         flipud(XX)];
    ciDD_back	= [ci_back(:,1); flipud(ci_back(:,2))];
end

%% - 5 - Finding the point of intersection and its uncertainty
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
    ImData(eb, kx, data);
    patch([-1e3, 1e3, 1e3, -1e3, -1e3],[KWIN(1), KWIN(1), KWIN(2), KWIN(2), KWIN(1)],...
        [0 1 0], 'linewidth', 1, 'facealpha', 0, 'edgecolor', [0 0 1]);
    img_props(); cbar_props(); colorbar off;
    axis([min(eb(:)), max(eb(:)),min(kx(:)), max(kx(:))]);
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
    ax.YAxisLocation = 'left';             % 'left' | 'right' | 'origin'
    % - Axis labels and limits
    xlabel('$$ \bf  E_{B} (eV) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
    % - Plotting the x- and y-axes
    axis([min(eb(:)), max(eb(:)), 0, 1.15]);
end


%% 6 - Appending data to MATLAB data structure
fitStr              = struct();
fitStr.vb_fit_args  = iparams;
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