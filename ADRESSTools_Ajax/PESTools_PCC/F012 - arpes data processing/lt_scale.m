% --- Function to perform a Scaling transformation
function [xDat, yDat] = lt_scale(xDat, yDat, scale_args)
% [xDat, yDat] = lt_scale(xDat, yDat, scale_args)
%   Function to perform a Scaling transformation on 2D data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xDat:           [N x M] array of the x-axis (theta/kx)
%   -   yDat:           [N x M] array of the y-axis (eb)
%   -   scale_args:     1x2 cell of {Sx, Sy}
%
%   OUT: 
%   -   xDat:           [N x M] array of the x-axis (theta/kx) after scaling
%   -   yDat:           [N x M] array of the y-axis (eb) after scaling

%% Default parameters
if nargin < 2; scale_args = {1, 1}; end
if isempty(scale_args); scale_args = {1, 1}; end

% disp('Scaling...')

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
xDat  = Sx .* (xDat - mean(xDat(:))) + mean(xDat(:));
yDat  = Sy .* (yDat - mean(yDat(:))) + mean(yDat(:));

end