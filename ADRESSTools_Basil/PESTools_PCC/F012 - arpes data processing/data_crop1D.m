function [xDat_crop, yDat_crop] = data_crop1D(xDat, yDat, xDat_lims)
% [xDat_crop, yDat_crop] = data_crop1D(xDat, yDat, xDat_lims)
%   This function crops the ARPES x-, y- and z-independent variables over
%   a given range. The function crops the most recently processed data and
%   ensures consistency across both the independent variables and data
%   matrix.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%
%   IN:
%   -   xDat:           [N x 1] array of the x-axis (theta/kx)
%   -   yDat:           [N x 1] array of the y-axis (eb)
% 	-   xDat_lims:      [1 x 2] row vector of x-axis limits
%
%   OUT:
%   -   xDat_crop:  	[A x B] array of the cropped x-axis (theta/kx)
%   -   yDat_crop:   	[A x B] array of the cropped y-axis (eb)

%% Default parameters
maxLim = 1e4;
if nargin < 3; xDat_lims=[-1 1]*maxLim; yDat_lims=[-1 1]*maxLim; end
if isempty(xDat_lims); xDat_lims=[-1 1]*maxLim;  end
% - Sorting the cropping limits in ascending order
xDat_lims = sort(xDat_lims);

%% 1 - Cropping the data
% - x-axis indices
[~, xIndxL]     = min(abs(xDat(:) - xDat_lims(1)));
[~, xIndxU]     = min(abs(xDat(:) - xDat_lims(2)));
x_indx          = [xIndxL xIndxU];
% - cropping x-axis
xDat_crop = xDat(x_indx(1):x_indx(2));
yDat_crop = yDat(x_indx(1):x_indx(2));

end
