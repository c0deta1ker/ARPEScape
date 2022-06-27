function loglog_props()
% loglog_props()
%   This function outlines a constant, global formatting for axes properties
%   and parameters for all figures plotted with the PESTools package. These
%   properties ensure that a log-log plot is made, with consistant
%   properties as 'gca_props()'.
%
%   REQ. FUNCTIONS: none
%
%   IN:     (none)
%
%   OUT:    (none)

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
ax.TickLength       = [0.02 0.025];
ax.XColor           = [0 0 0]; 
ax.YColor           = [0 0 0];
% Ruler properties
ax.XAxisLocation    = 'bottom';            % 'bottom' | 'top' | 'origin'
ax.YAxisLocation    = 'left';              % 'left' | 'right' | 'origin'
% Box Styling properties
ax.LineWidth        = 1.0;
ax.Box              = 'on'; 
ax.Layer            = 'Top';
% Making it a log-log plot
ax.XScale    = 'log';
ax.YScale    = 'log'; 
%% 2 - Plotting the x- and y-axes
xl = xlim; yl = ylim;
axis([xl(1), xl(2), yl(1), yl(2)]);
end