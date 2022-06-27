function fig = view_stat_hist(params, labels)
% fig = view_stat_hist(params, labels)
%   This function is used to plot the output paramaters of fitting functions, 
%   where you want to keep track of the value of a fit parameters vs run 
%   number. Here, a histogram is plotted of all the N solutions to see if
%   they exhibit a Gaussian / Uniform distribution, which may give some
%   information to the best solutions.
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
%
%   OUT:
%   -   fig:        MATLAB figure object with the histogram analysis plotted.

%% Default parameters
if nargin < 2; labels = []; end
if isempty(labels); labels = []; end
% -- Initialising the plot parameters
pp  = plot_props();

%% - 1 - INITIALISING VARIABLES
% -- Initialising the variables
P       = size(params, 2);      % number of independent parameters
N       = size(params, 1);      % number of independent runs
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'Statistical Analysis - Histogram of Fit Parameters');
fig.Position(3) = (P./2)*pp.fig5x4(1); 
fig.Position(4) = 0.7*pp.fig5x4(2);
%% - 2 - PLOTTING THE PARAMETERS VS N
for p = 1:P
    subplot(1,P,p); hold on;
    histogram(params(:,p), 'facecolor', pp.col.fit{p});
    % -- Formatting the figure
%     axis padded; 
    gca_props(); 
    if ~isempty(labels); title(labels{p}, 'fontweight', 'bold'); end
    ax = gca; ax.FontSize = 10;
end
end