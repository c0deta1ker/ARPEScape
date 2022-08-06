function img_props(colmap, type, plot_axes_lines)
% img_props(colmap, type)
%   This function outlines the consistent axes and figure
%   parameters to be used for 2D images. Each field can be edited to what
%   the user desires. Specific inputs are used to isolate the
%   different types of figures that are plotted. This allows for a constant
%   global formatting for 2D ARPES images within the PESTools package.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   colmap:     char of either: hot, parula, jet, gray, bone, hsv, cool, spring, 
%                   summer, autumn, winter, copper, pink, lines, colorcube, 
%                   prism, flag . 
%   -   type:       string of either "theta" or "kx" for the x-axis label
%
%   OUT: (none)

%% Default parameters
if nargin < 1; colmap = 'hot'; type = "kx"; end
if nargin < 2; type = "kx"; end
if nargin < 3; plot_axes_lines = 1; end
if isempty(colmap); colmap = 'hot'; end
if isempty(type);   type = "kx"; end
if isempty(plot_axes_lines);   plot_axes_lines = 1; end
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
ax.TickLength       = [0.01 0.025];
ax.XColor           = [0 0 0]; 
ax.YColor           = [0 0 0];
% Ruler properties
ax.XAxisLocation    = 'bottom';            % 'bottom' | 'top' | 'origin'
ax.YAxisLocation    = 'left';              % 'left' | 'right' | 'origin'
% Box Styling properties
ax.Color            = [0, 0, 0];
ax.LineWidth        = 1.2;
ax.Box              = 'on'; 
ax.Layer            = 'Top';
% Axis labels and limits
if type == "raw_tht" || type == "tht"
    xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex');
elseif type == "kx"
    xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
end
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
%% 2 - Defining the figure properties
% fig                 = gcf; 
% fig.Color           = [1 1 1]; 
% fig.InvertHardcopy  = 'off';
%% 3 - Colorbar properties
colormap(colmap);
%% 4 - Plotting the x- and y-axes
xl = xlim; yl = ylim;
axis([xl(1), xl(2), yl(1), yl(2)]);
% -- If true, plot axes lines and remove them from the legend
if plot_axes_lines == 1
    a = line([0 0], [-1e5, 1e5], 'Color', [1 1 1 0.5], 'LineWidth', 1., 'Linestyle', ':');
    b = line([-1e5, 1e5], [0 0], 'Color', [1 1 1 0.5], 'LineWidth', 1., 'Linestyle', ':');
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
end