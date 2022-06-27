function gca_props(plot_axes_lines)
% gca_props(plot_axes_lines)
%   This function outlines a constant, global formatting for axes properties
%   and parameters for all figures plotted with the PESTools package.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   plot_axes_lines:     if 1, will plot the x- and y-axes as black, dashed lines.
%
%   OUT:    (none)

%% Default parameters
if nargin < 1; plot_axes_lines = 1; end
if isempty(plot_axes_lines); plot_axes_lines = 1; end
%% 1 - Defining the axes properties
ax = gca;
% Font properties
ax.FontName         = 'Helvetica'; 
ax.FontWeight       = 'normal'; 
ax.FontSize         = 12;
% Tick properties
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick       = 'off'; 
ax.YMinorTick       = 'off';
ax.TickDir          = 'in';
ax.TickLength       = [0.015 0.025];
ax.XColor           = [0 0 0]; 
ax.YColor           = [0 0 0];
% Ruler properties
ax.XAxisLocation    = 'bottom';            % 'bottom' | 'top' | 'origin'
ax.YAxisLocation    = 'left';              % 'left' | 'right' | 'origin'
% Box Styling properties
ax.LineWidth        = 1.0;
ax.Box              = 'on'; 
ax.Layer            = 'Top';
%% 2 - Plotting the x- and y-axes
xl = xlim; yl = ylim;
axis([xl(1), xl(2), yl(1), yl(2)]);
% -- If true, plot axes lines and remove them from the legend
if plot_axes_lines == 1
    a = line([0 0], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    b = line([-1e5, 1e5], [0 0], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
end