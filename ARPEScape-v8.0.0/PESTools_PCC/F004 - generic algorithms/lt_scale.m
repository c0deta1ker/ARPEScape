function [xdat_lt, ydat_lt] = lt_scale(xdat, ydat, scale_args, plot_result)
% [xdat_lt, ydat_lt] = lt_scale(xdat, ydat, scale_args, plot_result)
%   Function to perform a Scaling transformation on 2D data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:           [nX×nY] array of the x-axis (theta/kx)
%   -   ydat:           [nX×nY] array of the y-axis (eb)
%   -   scale_args:     {1×2} cell of {Sx, Sy}
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT: 
%   -   xdat_lt:        [nX×nY] array of the x-axis (theta/kx) after scaling
%   -   ydat_lt:        [nX×nY] array of the y-axis (eb) after scaling

%% Default parameters
if nargin < 4; plot_result=0; end
if nargin < 3; scale_args = {1, 1}; end
if isempty(scale_args); scale_args = {1, 1}; end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
if length(scale_args) ~= 2; scale_args = {0, 0}; end

%% - 1 - Initialising the transformation parameters
% - Extracting the transformation values
Sx    = scale_args{1};
Sy    = scale_args{2};
% - Defining the transformation matrix
Smatrix = [...
    Sx, 0, 0;...
    0, Sy, 0;...
    0, 0, 1];
Stform  = affine2d(Smatrix);

%% 2 - Performing transformation operation
xdat_lt  = Sx .* (xdat - mean(xdat(:))) + mean(xdat(:));
ydat_lt  = Sy .* (ydat - mean(ydat(:))) + mean(ydat(:));

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