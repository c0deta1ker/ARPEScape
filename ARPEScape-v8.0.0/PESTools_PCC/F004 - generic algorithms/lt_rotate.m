function [xdat_lt, ydat_lt] = lt_rotate(xdat, ydat, rot_args, plot_result)
% [xdat_lt, ydat_lt] = lt_rotate(xdat, ydat, rot_args, plot_result)
%   Function to perform a Rotate transformation on 2D data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:           [nX×nY] array of the x-axis (theta/kx)
%   -   ydat:           [nX×nY] array of the y-axis (eb)
%   -   rot_args:       {1×1} cell of {Rtheta}
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT: 
%   -   xdat_lt:        [nX×nY] array of the x-axis (theta/kx) after rotation
%   -   ydat_lt:        [nX×nY] array of the y-axis (eb) after rotation

%% Default parameters
if nargin < 4; plot_result=0; end
if nargin < 3; rot_args = {0}; end
if isempty(rot_args); rot_args = {0}; end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
if length(rot_args) ~= 1; rot_args = {0}; end

%% - 1 - Initialising the transformation parameters
% - Extracting the transformation values
Rtheta    = rot_args{1};
% - Defining the 2D rotation matrix
R3Dmatrix = [...
    cos(deg2rad(Rtheta)),  sin(deg2rad(Rtheta));...
    -sin(deg2rad(Rtheta)), cos(deg2rad(Rtheta))];
%% 2 - Performing transformation operation
if length(xdat(:)) == length(ydat(:))
    XY      = [xdat(:) ydat(:)];  	% Create Matrix Of Vectors
    rotXY   = XY*R3Dmatrix';        % Multiply the vectors by the rotation matrix
    xdat_lt = reshape(rotXY(:,1), size(xdat,1), []);
    ydat_lt = reshape(rotXY(:,2), size(ydat,1), []);
else
    xdat = repmat(xdat, [length(ydat), 1]);
    ydat = repmat(ydat, [1, length(xdat)]);
    xdat_lt   = xdat*R3Dmatrix;
    ydat_lt   = ydat*R3Dmatrix;
end

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