function cb = cbar_props(colmap, varargin)
% cb = cbar_props(varargin)
%   This function outlines the consistent parameters to be used for the
%   colorbar settings. The user can control the position and side of the
%   tick marks, otherwise, the settings below are used for all figures
%   plotted with the PESTools package. This allow for a constant, global 
%   formatting for colorbar properties.
%
%   IN:
%   -   colmap:     char of either: hot, parula, jet, gray, bone, hsv, cool, spring, 
%                   summer, autumn, winter, copper, pink, lines, colorcube, prism, flag . 
%   -   varargin:   ColProps arguments: ColorBar properties control the appearance and behavior of a Colorbar object, specified as a string scalar or character vector.
%
%   OUT:
%   -   cb:         ColorBar properties.

%% Default parameters
if nargin < 1; colmap = 'hot'; end
if nargin < 2; varargin = {}; end
if isempty(colmap); colmap = 'hot'; end
if isempty(varargin); varargin = {}; end
%% 1 - Defining the colormap
colormap(colmap);
%% 2 - Defining the colorbar
cb = colorbar;
% Colorbar font properties
cb.FontName             = 'Segoe UI';  
cb.FontWeight           = 'normal'; 
cb.FontSize             = 8;
% Colorbar tick properties
cb.TickLength           = 0.02; 
cb.TickDirection        = 'both';
cb.YAxisLocation        = 'right';
% Colorbar position properties
cb.Location             = 'eastoutside';  % 'north' | 'south' | 'east' | 'west' | 'northoutside' | 'southoutside' | 'eastoutside' | 'westoutside'
cb.Position             = [0.925 0.125 0.025 0.175];
% Colorbar box properties
cb.Color                = [0 0 0]; 
cb.Box                  = 'on'; 
cb.LineWidth            = 0.75;
cb.Ticks                = sort([cb.Limits(1), 0.5*(cb.Limits(2)), cb.Limits(2)]);
%% 3 - Setting the colorbar properties as needed
if ~isempty(varargin); set(cb, varargin{:}); end
end