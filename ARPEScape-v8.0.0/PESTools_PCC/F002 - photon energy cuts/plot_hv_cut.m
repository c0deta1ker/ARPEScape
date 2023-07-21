function h = plot_hv_cut(hvCut, plot_args)
% h = plot_hv_cut(hvCut, plot_args)
%   This is a function that plots the planar Brilluoin zone 
%   extracted from the 'get_hv_cut()' function.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   hvCut:          MATLAB data structure containing hc cut data from 'get_hv_cut()'
%   -   plot_args:      plot arguments for the Brilluoin zone slice
%
%   OUT: (none)

%% Default parameters
if nargin < 2; plot_args = []; end
if isempty(plot_args); plot_args = []; end
%% - 1 - Initialising the transformation parameters
hold on;
ax = gca;
XLim = ax.XLim; YLim = ax.YLim; 
col = lines(2);
% - For a single photon energy
if length(hvCut.hv) == 1
    % -- Plotting the ARPES scan line
    plot(hvCut.kx, hvCut.kz, '-', 'linewidth', 2, 'color', col(1,:));
    % -- Adding text
    text(max(hvCut.kx(:)), min(hvCut.kz(:)), sprintf("%.0f eV", hvCut.hv),...
        'color', col(1,:), 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
% - For a range of photon energies
elseif length(hvCut.hv) == 2
    % -- Plotting the ARPES scan plane which includes the hvCut.kz broadening
    patch([hvCut.kx(:,2); flipud(hvCut.kx(:,1))], [hvCut.kz(:,2); flipud(hvCut.kz(:,1))], col(2,:), 'edgecolor', 'none', 'facealpha', 0.5);
    % -- Plotting the ARPES scan line
    plot(hvCut.kx(:,1), hvCut.kz(:,1), '-', 'linewidth', 2, 'color', col(1,:));
    plot(hvCut.kx(:,2), hvCut.kz(:,2), '-', 'linewidth', 2, 'color', col(2,:));
    % -- Adding text
    text(max(hvCut.kx(:,1)), min(hvCut.kz(:,1)), sprintf("%.0f eV", hvCut.hv(1)),...
        'color', col(1,:), 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
    text(max(hvCut.kx(:,2)), min(hvCut.kz(:,2)), sprintf("%.0f eV", hvCut.hv(2)),...
        'color', col(2,:), 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
end
axis([XLim, YLim]);
end