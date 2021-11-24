function pp = plot_props()
% pp = plot_props()
%   This function outlines the consistent colors, markers and lines to be
%   used when plotting the figures associated with all the state fitting
%   analysis. Whenever a figure is plotted using the PESTools package, you
%   will most likely see this function being called. 
%
%   REQ. FUNCTIONS: none
%
%   IN:     (none)
%
%   OUT:
%   -   pp:     MATLAB structure that contains all of the figure and plot properties.

% - Initialising plot property structure
pp = struct();

%% FIGURE SIZE PROPERTIES
pp.fontsize     = 12;
% -- LANDSCAPE FIGURE SIZE
pp.fig3x2       = [600, 400];
pp.fig4x3       = [800, 600];
pp.fig4x5       = [400, 500];   % best size in my opinion for portrait
pp.fig16x9  	= [640, 360];
pp.fig16x10 	= [1280, 800];
pp.fig16x15 	= [1280, 1200];
% -- PORTRAIT FIGURE SIZE
pp.fig2x3       = [400, 600];
pp.fig3x4       = [600, 800];
pp.fig5x4       = [500, 400];   % best size in my opinion for landscape
pp.fig9x16  	= [360, 640];
pp.fig10x16 	= [800, 1280];
pp.fig15x16 	= [1200, 1280];
pp.fig13x12   	= [600, 650];
% -- SQUARE FIGURE SIZE
pp.fig4x4       = [400 400];
pp.fig6x6       = [600 600];
pp.fig8x8       = [800 800];
% -- FIGURE POSITION
pp.figpos       = [100 100];

%% GENERAL PROPERTIES
pp.lwidth       = 1.0;
pp.llwidth      = 2.0;
pp.msize        = 9.0;
pp.mmsize       = 12.0;
pp.falpha       = 0.50;
pp.ffalpha      = 0.80;
pp.csize        = 3.0;
pp.ccsize       = 6.0;
pp.bzlwidth     = 2.0;
pp.clims        = [0, 1];

%% DEFINING THE COLORS SCHEME
% - FIT COMPONENT COLORS
pp.col.dat = {...
    [0.9, 0.1, 0.1],...                 %% RED FOR EXPERIMENTAL DATA
    [0, 0, 0]...                        %% BLACK FOR BEST FIT
    };
% -- CUT COMPONENT COLOURS
pp.col.cut      = {...
    [0.0, 0.0, 0.0],...                 %% BLACK    (FDD)
    [0.2, 0.6, 0.2],...                 %% GREEN    (EDC)
    [0.2, 0.2, 0.8],...                 %% BLUE     (MDC)
    [0.1, 0.1, 0.6],...                 %% BLUE     (ROI)
    };
% - FIT COMPONENT COLORS
pp.col.fit = num2cell([0.25,0.25,0.25;lines(20)], 2);

%% DEFINING THE MARKER STYLES
% -- CUT COMPONENT COLOURS
pp.mstyle.cut      = {...
    'x',...    % (FDD)
    's',...    % (EDC)
    '^',...    % (MDC)
    'x',...    % (ROI)
    };
% -- FIT COMPONENT COLOURS
pp.mstyle.fits     = {...
    'x',...    % (n = 1)
    '+',...    % (n = 2)
    'o',...    % (n = 3)
    '^',...    % (n = 4)
    '.',...    % (n = 5)
    's',...    % (n = 6)
    'p',...    % (n = 7)
    };

%% MISCELLANEOUS PLOT PROPERTIES
% -- High-symmtry ticks for lattice vectors in silicon
pp.si_ebkticks = {...
    sort(-116:1.16:116),...                         % GX
    sort(-116:0.5*1.16:116),...                     % XW
    sort([-5*1.64:1.64:5*1.64, 0.41, -0.41]),...  	% GKX
    };
end