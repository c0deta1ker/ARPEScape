function ax = img_props(plot_axes_lines, varargin)
% ax = img_props(plot_axes_lines, varargin)
%   This function outlines a constant, global formatting for axes properties
%   and parameters for all 2D images plotted with the PESTools package. In
%   particular, this is used for all ARPES data images.
%
%   IN:
%   -   plot_axes_lines:    if 1, will plot the x- and y-axes as black, dashed lines.
%   -   varargin:           AxProps arguments: Axes properties control the appearance and behavior of an Axes object, specified as a string scalar or character vector.
%
%   OUT:
%   -   ax:                 Axes properties.

%% Default parameters
if nargin < 1; plot_axes_lines = 1; end
if nargin < 2; varargin = {}; end
if isempty(plot_axes_lines); plot_axes_lines = 1; end
if isempty(varargin); varargin = {}; end
%% 1 - Defining the axes properties
ax = gca;
% Font properties
ax.FontName         = 'Segoe UI'; 
ax.FontWeight       = 'normal'; 
ax.FontSize         = 10;
% Tick properties
ax.XMinorTick       = 'off'; 
ax.YMinorTick       = 'off';
ax.TickDir          = 'both';
ax.XColor           = [0 0 0]; 
ax.YColor           = [0 0 0];
% Ruler properties
ax.XAxisLocation    = 'bottom';             % 'bottom' | 'top' | 'origin'
ax.YAxisLocation    = 'left';               % 'left' | 'right' | 'origin'
% Box Styling properties
ax.Color            = [0 0 0];
ax.Box              = 'off';                % 'on' | 'off'
ax.LineWidth        = 0.75; 
ax.Layer            = 'Top';
% Axis Scale
ax.XScale           = 'linear';             % 'linear' | 'log'
ax.YScale           = 'linear';             % 'linear' | 'log' 
%% 2 - Plotting the x- and y-axes
xl = xlim; yl = ylim;
axis([xl(1), xl(2), yl(1), yl(2)]);
% -- If true, plot axes lines and remove them from the legend
if plot_axes_lines == 1
    a = xline(0, 'Color', [1 1 1], 'LineWidth', 1.00, 'Linestyle', ':');
    b = yline(0, 'Color', [1 1 1], 'LineWidth', 1.00, 'Linestyle', ':');
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
%% 3 - Setting the axes properties as needed
if ~isempty(varargin); set(ax, varargin{:}); end
end