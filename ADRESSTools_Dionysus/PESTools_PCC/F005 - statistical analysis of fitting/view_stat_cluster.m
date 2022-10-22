function fig = view_stat_cluster(params, labels, n_order)
% fig = view_stat_cluster(params, labels, n_order)
%   This function is used to plot the output paramaters of fitting functions, 
%   where you want to keep track of the value of a fit parameters vs run 
%   number. This can be used to check for any large deviations of the fit
%   parameters and whether you should impose any constraints if some
%   solutions are unphysical. The cluster analysis plots consist of 3
%   layers; (1) the first is a pair-wise scatter plot of the fit
%   parameters, allowing you to observe how repeatable and convergent each
%   solution is; (2) the second is a colormap which is overlaid on the
%   background of each pair-wise plot, showing you the frequency of each
%   fit parameter. This allows you to identify the most probable solution.
%   (3) if 'n_order' does not equal zero, it will plot a 2D Gaussian
%   function around the pair-wise clusters to give an estimate of their
%   mean position and standard deviation.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   params:     (N x P) array of P fit parameters converged with 
%                   N independent runs. For a system in which some
%                   parameters are common, you can define the params as a 
%                   {1 x P} cell array, each of which has an (N x C) array
%                   of C common fit parameters with N independent runs.
%   -   labels:     (1 x P) cell array that contains the labels of each fit
%                   parameter.
%   -   n_order:    scalar of the number of Gaussians to use when trying to
%                   fit the cluster analysis.
%
%   OUT:
%   -   fig:        MATLAB figure object with the cluster analysis plotted.

%% Default parameters
if nargin < 3; n_order = 0; end
if nargin < 2; labels = []; end
if isempty(labels); labels = []; end
if isempty(n_order); n_order = 0; end
% -- Initialising the plot parameters
pp  = plot_props();

%% - 1 - INITIALISING VARIABLES
% -- Initialising the variables
P       = size(params, 2);      % number of independent parameters
N       = size(params, 1);      % number of independent runs

%% - 2 - PLOTTING THE BEST FIT PARAMETERS AS A FIGURE MATRIX
for p = 1:P
    % -- Initialising the figure
    fig(p) = figure(); set(fig(p), 'Name', 'Statistical Analysis - Parameters vs Run Number');
    fig(p).Position(3) = P*pp.fig4x4(1); 
    fig(p).Position(4) = 1.0*pp.fig4x4(2);
    for j = 1:P
        % -- Plotting the best fit parameters as a cluster plot
        subplot(1,P,j); hold on;
        plot(params(:,p), params(:,j), 'k.', 'color', pp.col.fit{p}, 'markerfacecolor', pp.col.fit{p}, 'markersize', 10);
        % -- Plotting the best fit parameters as a coloured histogram
        N = hist3([params(:,p), params(:,j)], 'CDataMode','auto','FaceColor','interp');
        N_pcolor = N'; N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
        xl = linspace(min(params(:,p)),max(params(:,p)),size(N_pcolor,2));  % Columns of N_pcolor
        yl = linspace(min(params(:,j)),max(params(:,j)),size(N_pcolor,1));  % Rows of N_pcolor
        h = pcolor(xl,yl,N_pcolor);
        h.ZData = -1*ones(size(N_pcolor));      % Shifting the map to be below data points
        % -- Plotting the best fit Gaussian around the clusters
        if n_order ~= 0
            if p ~= j
                A           = params(:,[p,j]);
                GMModel     = fitgmdist(A,n_order);
                gmPDF       = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
                g = gca; fcontour(gmPDF,[g.XLim g.YLim]);
            end
        end
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