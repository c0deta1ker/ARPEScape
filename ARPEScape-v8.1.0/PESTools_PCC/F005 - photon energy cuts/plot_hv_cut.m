function h = plot_hv_cut(hvCuts, varargin)
% h = plot_hv_cut(hvCut, plot_args)
%   This is a function that plots the planar Brilluoin zone 
%   extracted from the 'get_hv_cut()' function.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   hvCuts:         1xN cell array of the MATLAB data structure containing hv cut data from 'get_hv_cut()'
%   -   plot_args:      plot arguments for the Brilluoin zone slice
%
%   OUT: (none)

%% Default parameters
if nargin < 2; varargin = {}; end
if isempty(varargin); varargin = {}; end
%% - 1 - Plotting the photon energy cuts
hold on;
ax = gca;
XLim = ax.XLim; YLim = ax.YLim;
axis([XLim, YLim]);
% -- Plotting grid lines of photon energies
hvGrid = 200:200:3000;
for i = 1:length(hvGrid)
    hvStruct = get_hv_cut(hvGrid(i));
    plot(hvStruct.kx, hvStruct.kz, 'r:', 'linewidth', 1);
    if hvStruct.hv == 200 || hvStruct.hv == 500 || hvStruct.hv == 1000 || hvStruct.hv == 2000
        text(max(hvStruct.kx(:)), min(hvStruct.kz(:)), sprintf("%.0f eV", hvStruct.hv), 'color', 'k', 'Fontsize', 8, 'horizontalalignment', 'left');
    end
end
% -- Plotting the required cuts
if isstruct(hvCuts)
    col = [0 0 0];
    plot(hvCuts.kx, hvCuts.kz, '-', 'color', col, 'linewidth', 2, varargin{:});
    text(max(hvCuts.kx(:)), min(hvCuts.kz(:)), sprintf("%.0f eV", hvCuts.hv), 'color', col, 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
else
    N = length(hvCuts); col = lines(N);
    for i = 1:N
        plot(hvCuts{i}.kx, hvCuts{i}.kz, '-', 'color', col(i,:), 'linewidth', 2, varargin{:});
        text(max(hvCuts{i}.kx(:)), min(hvCuts{i}.kz(:)), sprintf("%.0f eV", hvCuts{i}.hv), 'color', col(i,:), 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
    end
end
end