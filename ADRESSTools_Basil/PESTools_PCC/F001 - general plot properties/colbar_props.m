function colbar_props(cbar_pos, cbar_tickside)
% colbar_props(cbar_pos, cbar_tickside)
%   This function outlines the consistent parameters to be used for the
%   colorbar settings. The user can control the position and side of the
%   tick marks, otherwise, the settings below are used for all figures
%   plotted with the PESTools package. This allow for a constant, global 
%   formatting for colorbar properties.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   cbar_pos:           1 x 4 vector of colorbar position in normalised units [left, bottom, width, height]
%   -   cbar_tickside:   	char of either 'left' or 'right' for the y-axis ticks
%
%   OUT: (none)

%% Default parameters
if nargin < 1; cbar_pos = [0.925 0.125 0.025 0.175]; cbar_tickside='right'; end
if nargin < 2; cbar_tickside = 'right'; end
if isempty(cbar_pos); cbar_pos = [0.925 0.125 0.025 0.175]; end
if isempty(cbar_tickside); cbar_tickside='right'; end
%% 1 - Defining the colorbar
cb = colorbar;
% Colorbar font properties
cb.FontName             = 'Helvetica'; 
cb.FontSize             = 9;
% Colorbar tick properties
cb.TickLabelInterpreter = 'latex';
cb.TickLength           = 0.04; 
cb.TickDirection        = 'both';
cb.YAxisLocation        = cbar_tickside;
% Colorbar position properties
cb.Position             = cbar_pos;
% Colorbar box properties
cb.Color                = [0 0 0]; 
cb.Box                  = 'on'; 
cb.LineWidth            = 1.0;
if cb.Limits(1) < 0 || cb.Limits(1) == 0
        cb.Ticks = sort([0.98*cb.Limits(1), 0.5*(cb.Limits(2)), 0.98*cb.Limits(2)]);
else;   cb.Ticks = sort([1.02*cb.Limits(1), 0.5*(cb.Limits(2)), 0.98*cb.Limits(2)]); 
end
end