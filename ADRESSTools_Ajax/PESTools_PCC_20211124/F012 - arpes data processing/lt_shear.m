% --- Function to perform a Shear transformation
function [xDat, yDat] = lt_shear(xDat, yDat, shear_args)
% [xDat, yDat] = lt_shear(xDat, yDat, shear_args)
%   Function to perform a Shear transformation on 2D data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xDat:           [N x M] array of the x-axis (theta/kx)
%   -   yDat:           [N x M] array of the y-axis (eb)
%   -   shear_args:     1x2 cell of {Sx, Sy}
%
%   OUT: 
%   -   xDat:           [N x M] array of the x-axis (theta/kx) after shearing
%   -   yDat:           [N x M] array of the y-axis (eb) after shearing

%% Default parameters
if nargin < 2; shear_args = {0, 0}; end
if isempty(shear_args); shear_args = {0, 0}; end

% disp('Shearing...')

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
if      Shx == 0; x_shifts = zeros(size(yDat));
else;   x_shifts = shear_fnc(Shx, 0, yDat);
end
% - Applying y-shear corrections
if      Shy == 0; y_shifts = zeros(size(xDat));
else;   y_shifts = shear_fnc(Shy, 0, xDat);
end
% - Extracting the initial mid-point of the data
Ox = mean(xDat(:));
Oy = mean(yDat(:));
% - Execute shear transformation
xDat = (xDat + x_shifts);
yDat = (yDat + y_shifts);
% - Recentering the data after shearing
xDat = xDat + (Ox - mean(xDat(:)));
yDat = yDat + (Oy - mean(yDat(:)));

end