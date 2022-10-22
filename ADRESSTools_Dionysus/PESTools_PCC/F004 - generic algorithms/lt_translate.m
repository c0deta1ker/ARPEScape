function [xDat, yDat] = lt_translate(xDat, yDat, trans_args)
% [xDat, yDat] = lt_translate(xDat, yDat, trans_args)
%   Function to perform a Translation transformation on 2D data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xDat:           [N×M] array of the x-axis (theta/kx)
%   -   yDat:           [N×M] array of the y-axis (eb)
%   -   trans_args:     {1×2} cell of {Tx, Ty}
%
%   OUT: 
%   -   xDat:           [N×M] array of the x-axis (theta/kx) after translation
%   -   yDat:           [N×M] array of the y-axis (eb) after translation

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
xDat  = xDat - Tx;
yDat  = yDat - Ty;

end