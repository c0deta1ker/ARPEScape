function [xdat_lt, ydat_lt] = lt_translate(xdat, ydat, trans_args, plot_result)
% [xdat_lt, ydat_lt] = lt_translate(xdat, ydat, trans_args, plot_result)
%   Function to perform a Translation transformation on 2D data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:           [nX×nY] array of the x-axis (theta/kx)
%   -   ydat:           [nX×nY] array of the y-axis (eb)
%   -   trans_args:     {1×2} cell of {Tx, Ty}
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT: 
%   -   xdat_lt:        [nX×nY] array of the x-axis (theta/kx) after translation
%   -   ydat_lt:        [nX×nY] array of the y-axis (eb) after translation

%% Default parameters
if nargin < 3; trans_args = {0, 0}; end
if isempty(trans_args); trans_args = {0, 0}; end
%% Validity checks on the input parameters
if length(trans_args) ~= 2; trans_args = {0, 0}; end

%% - 1 - Initialising the transformation parameters
% - Extracting the transformation values
Tx    = trans_args{1};
Ty    = trans_args{2};
% - Defining the transformation matrix
Tmatrix = [...
    1, 0, 0;...
    0, 1, 0;...
    Tx, Ty, 1];
Ttform  = affine2d(Tmatrix);

%% 2 - Performing transformation operation
xdat_lt  = xdat + Tx;
ydat_lt  = ydat + Ty;

%% -- For Debugging
if plot_result == 1
    fig = figure(); 
    fig.Position(3) = 400*2; 
    fig.Position(4) = 400*0.9; 
    subplot(121); hold on;
    h1 = pcolor(xdat, ydat, ydat.*0); set(h1,'EdgeColor','None','FaceColor','Flat');
    img_props();
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
    subplot(122); hold on;
    h2 = pcolor(xdat_lt, ydat_lt, ydat_lt.*0); set(h2,'EdgeColor','None','FaceColor','Flat');
    img_props();
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    axis([min(xdat_lt(:)), max(xdat_lt(:)), min(ydat_lt(:)), max(ydat_lt(:))]);
end
end