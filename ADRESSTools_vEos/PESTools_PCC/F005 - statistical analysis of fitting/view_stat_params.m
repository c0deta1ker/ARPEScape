function fig = view_stat_params(params, labels)
% fig = view_stat_params(params, labels)
%   This function is used to plot the output parameters of fitting functions, 
%   where you want to keep track of the value of a fit parameters vs run 
%   number. This can be used to check for any large deviations of the fit
%   parameters and whether you should impose any constraints if some
%   solutions are unphysical.
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
%   -   fig:        MATLAB figure object with the fit parameter analysis plotted.

%% Default parameters
if nargin < 2; labels = []; end
if isempty(labels); labels = []; end
% -- Initialising the plot parameters
pp  = plot_props();

%% - 1 - INITIALISING VARIABLES
if iscell(params)
    P  	= length(params);       % number of independent parameters
    N  	= size(params{1}, 1); 	% number of independent runs
else
    P  	= size(params, 2);      % number of independent parameters
    N  	= size(params, 1);      % number of independent runs
end
nruns   = 1:N;                  % number of independent runs as a vector
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'Statistical Analysis - Parameters vs Run Number');
fig.Position(3) = (P./2)*pp.fig5x4(1); 
fig.Position(4) = 1.25*pp.fig5x4(2);

%% - 2 - PLOTTING THE FIT PARAMETERS VS NUMBER OF RUNS
for p = 1:P
    subplot(1,P,p); hold on;
    % -- Plotting the outcome for the XPS fitting 
    if iscell(params)
        C  	= size(params{p}, 2); 	% number of independent common fit parameters.
        for c = 1:C
            plot(params{p}(:,c), nruns, 'k.-', 'color', pp.col.fit{c}, 'markerfacecolor', pp.col.fit{c}, 'linewidth', 1.25);
            errorbar(mean(params{p}(:,c)), mean(nruns(:)), 0, 0, 3*std(params{p}(:,c)), 3*std(params{p}(:,c)), 'ks-',...
                'color', [0 0 0], 'markerfacecolor', pp.col.fit{c}, 'markersize', 12, 'linewidth', 1.5);
        end
    % -- Plotting the outcome for the ARPES fitting
    else
        plot(params(:,p), nruns, 'k.-', 'color', pp.col.fit{p}, 'markerfacecolor', pp.col.fit{p}, 'linewidth', 1.25);
        errorbar(mean(params(:,p)), mean(nruns(:)), 0, 0, 3*std(params(:,p)), 3*std(params(:,p)), 'ks-',...
            'color', [0 0 0], 'markerfacecolor', pp.col.fit{p}, 'markersize', 12, 'linewidth', 1.5);
    end
    % -- Formatting the figure
%     axis padded;
    gca_props(); 
    if ~isempty(labels); title(labels{p}, 'fontweight', 'bold'); end
    ax = gca; ax.FontSize = 10;
end
end