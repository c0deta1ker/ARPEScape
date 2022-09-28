function fig = view_pes_data(pesStr)
% fig = view_pes_data(pesStr)
%   This function plots the XPS data in the form of I(X). The
%   input can either be a single XPS data structure, or a cell-array of XPS
%   structures, each of which can be plotted and overlaid on top of one
%   another. Use this as the general function to use to view all XPS data
%   during the data processing / analysis stage.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   pesStr:     	data structure (or cell array of structures) of the XPS data.
%
%   OUT:
%   -   fig:           	MATLAB figure object with the ARPES data plotted.

%% - 1 - Plotting the XPS figure
% -- Initialising the plot properties
pp  = plot_props();
% -- Initialising the figure
fig = figure(); set(fig, 'Name', 'XPS Data from ADRESS');
fig.Position(3) = 1.25*pp.fig5x4(1); 
fig.Position(4) = pp.fig5x4(2);
hold on;
% -- Filing through all potential XPS structures in a cell array
lgnd = {};
if length(pesStr) == 1 && ~iscell(pesStr)
    plot(pesStr.xdat, pesStr.ydat, 'k-', 'color', [0 0 1], 'linewidth', pp.llwidth);
    axis([min(pesStr.xdat(:)), max(pesStr.xdat(:)),0, 1.15*max(pesStr.ydat(:))]);
    lgnd{end+1} = string(pesStr.hv) + " eV";
    if isfield(pesStr, 'FileName'); title(sprintf(string(pesStr.FileName)), 'interpreter', 'none', 'fontsize', 9); end
else
    cols = jet(length(pesStr)+2);
    for i = 1:length(pesStr)
        plot(pesStr{i}.xdat, pesStr{i}.ydat, 'k-', 'color', cols(i,:), 'linewidth', pp.llwidth);
        lgnd{end+1} = string(pesStr{i}.hv) + " eV";
        max_ydat(i) = max(pesStr{i}.ydat(:));
    end
    axis([min(pesStr{1}.xdat(:)), max(pesStr{1}.xdat(:)), 0, 1.15*max(max_ydat)]);
    if isfield(pesStr, 'FileName'); title(sprintf(string(pesStr{1}.FileName)), 'interpreter', 'none', 'fontsize', 9); end
end
% - Formatting the figure
gca_props(); grid on;
% legend(lgnd,'Location','best');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');

end