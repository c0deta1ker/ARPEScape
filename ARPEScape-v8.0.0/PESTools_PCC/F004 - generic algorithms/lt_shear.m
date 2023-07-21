function [xdat_lt, ydat_lt] = lt_shear(xdat, ydat, shear_args, plot_result)
% [xdat_lt, ydat_lt] = lt_shear(xdat, ydat, shear_args, plot_result)
%   Function to perform a Shear transformation on 2D data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:           [nX×nY] array of the x-axis (theta/kx)
%   -   ydat:           [nX×nY] array of the y-axis (eb)
%   -   shear_args:     {1×2} cell of {Shx, Shy}
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT: 
%   -   xdat_lt:        [nX×nY] array of the x-axis (theta/kx) after shearing
%   -   ydat_lt:        [nX×nY] array of the y-axis (eb) after shearing

%% Default parameters
if nargin < 4; plot_result=0; end
if nargin < 3; shear_args = {0, 0}; end
if isempty(shear_args); shear_args = {0, 0}; end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
if length(shear_args) ~= 2; shear_args = {0, 0}; end

%% - 1 - Initialising the transformation parameters
% - Extracting the transformation values
Shx    = shear_args{1};
Shy    = shear_args{2};
% - Defining the transformation matrix
Smatrix = [...
    1, Shy, 0;...
    Shx, 1, 0;...
    0, 0, 1];
Stform  = affine2d(Smatrix);

%% 2 - Performing transformation operation
shear_fnc = @ (slope, intercept, y) (y - intercept) ./ slope;
% - Applying x-shear corrections
if      Shx == 0; x_shifts = zeros(size(ydat));
else;   x_shifts = shear_fnc(Shx, 0, ydat);
end
% - Applying y-shear corrections
if      Shy == 0; y_shifts = zeros(size(xdat));
else;   y_shifts = shear_fnc(Shy, 0, xdat);
end
% - Extracting the initial mid-point of the data
Ox = mean(xdat(:));
Oy = mean(ydat(:));
% - Execute shear transformation
xdat_lt = (xdat + x_shifts);
ydat_lt = (ydat + y_shifts);
% - Recentering the data after shearing
xdat_lt = xdat_lt + (Ox - mean(xdat_lt(:)));
ydat_lt = ydat_lt + (Oy - mean(ydat_lt(:)));

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