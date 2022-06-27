function fig = view_stat_cluster_prctile(params, labels, chisq, prc)
% fig = view_stat_cluster_prctile(params, labels, chisq, prc)
%   This function plots the same figure as in view_stat_cluster(), however,
%   it isolates the clusters into two segments; either above or below the
%   percentile [prc] value of the chi-squared [chisq] matrix. This allows
%   you to identify the top 10% of the optimal fits, allowing you to select
%   the best fit parameters more precisely.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   params:     (N x P) array of P fit parameters converged with 
%                   N independent runs. For a system in which some
%                   parameters are common, you can define the params as a 
%                   {1 x P} cell array, each of which has an (N x C) array
%                   of C common fit parameters with N independent runs.
%   -   labels:     (1 x P) cell array that contains the labels of each fit parameter.
%   -   chisq:      (1 x M) vector that contains the chi-squared values for each one of the converged fits.
%   -   prc:        scalar value between [0, 100] that thresholds all the fit 
%                   parameters within a percentile, 0 < prc, of the chi-
%                   squared values. This allows you to isolate the best
%                   solutions for diagnostic testing.
%
%   OUT:
%   -   fig:        MATLAB figure object with the cluster analysis plotted.

%% Default parameters
if nargin < 4; prc = 10; end
if isempty(chisq); error('Not enough arguments; make sure you input chisq.'); end
if isempty(labels); labels = []; end
if isempty(prc); prc = 10; end
% -- Initialising the plot parameters
pp  = plot_props();

%% - 1 - INITIALISING VARIABLES
% -- Initialising the variables
P       = size(params, 2);      % number of independent parameters
% -- Extracting the percentile value 
chisq_prctile   = prctile(chisq, prc);
prctile_indx    = chisq < chisq_prctile;
raw_indx        = ~prctile_indx;

%% - 2 - PLOTTING THE BEST FIT PARAMETERS AS A FIGURE MATRIX
for p = 1:P
    % -- Initialising the figure
    fig(p) = figure();
    fig(p).Position(3) = P*pp.fig4x4(1); 
    fig(p).Position(4) = 1.0*pp.fig4x4(2);
    for j = 1:P
        % -- Plotting the best fit parameters as a cluster plot
        subplot(1,P,j); hold on;
        plot(params(raw_indx,p), params(raw_indx,j), 'r.', 'markersize', 6);
        plot(params(prctile_indx,p), params(prctile_indx,j), 'g.', 'markersize', 10);
        % -- Plotting the best fit parameters as a coloured histogram
        N = hist3([params(:,p), params(:,j)], 'CDataMode','auto','FaceColor','interp');
        N_pcolor = N'; N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
        xl = linspace(min(params(:,p)),max(params(:,p)),size(N_pcolor,2));  % Columns of N_pcolor
        yl = linspace(min(params(:,j)),max(params(:,j)),size(N_pcolor,1));  % Rows of N_pcolor
        h = pcolor(xl,yl,N_pcolor);
        h.ZData = -1*ones(size(N_pcolor));      % Shifting the map to be below data points
        % -- Formatting the figure
        xlabel(labels{p}, 'interpreter', 'latex');
        ylabel(labels{j}, 'interpreter', 'latex');
        if min(params(:,p)) ~= max(params(:,p)) && min(params(:,j)) ~= max(params(:,j))
            axis([min(params(:,p)),max(params(:,p)),min(params(:,j)),max(params(:,j))]);
        end
        gca_props(); ax = gca; ax.FontSize = 9; axis square;
    end
end
end