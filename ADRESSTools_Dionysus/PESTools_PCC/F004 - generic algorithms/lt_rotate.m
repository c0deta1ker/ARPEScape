function [xDat, yDat] = lt_rotate(xDat, yDat, rot_args)
% [xDat, yDat] = lt_rotate(xDat, yDat, rot_args)
%   Function to perform a Rotate transformation on 2D data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xDat:           [N×M] array of the x-axis (theta/kx)
%   -   yDat:           [N×M] array of the y-axis (eb)
%   -   rot_args:       {1×1} cell of {Rtheta}
%
%   OUT: 
%   -   xDat:           [N×M] array of the x-axis (theta/kx) after rotation
%   -   yDat:           [N×M] array of the y-axis (eb) after rotation

%% Default parameters
if nargin < 3; rot_args = {0}; end
if isempty(rot_args); rot_args = {0}; end
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
XY      = [xDat(:) yDat(:)];  	% Create Matrix Of Vectors
rotXY   = XY*R3Dmatrix';        % Multiply the vectors by the rotation matrix
xDat = reshape(rotXY(:,1), size(xDat,1), []);
yDat = reshape(rotXY(:,2), size(yDat,1), []);

end