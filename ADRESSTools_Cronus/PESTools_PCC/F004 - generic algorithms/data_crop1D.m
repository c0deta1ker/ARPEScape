function [xdat_crop, ydat_crop] = data_crop1D(xdat, ydat, xdat_lims)
% [xdat_crop, ydat_crop] = data_crop1D(xdat, ydat, xdat_lims)
%   This function crops 1D data along the x-axis (domain) and y-axis
%   (range). The user defines the limits along the x-axis. 
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   xdat:           [N×1] array of the x-axis
%   -   ydat:           [N×1] array of the y-axis
% 	-   xdat_lims:      [1×2] row vector of x-axis limits
%
%   OUT:
%   -   xdat_crop:  	[A×B] array of the cropped x-axis
%   -   ydat_crop:   	[A×B] array of the cropped y-axis

%% Default parameters
if nargin < 3; xdat_lims=[min(xdat(:)), max(xdat(:))]; end
if isempty(xdat_lims); xdat_lims=[min(xdat(:)), max(xdat(:))];  end
% - Sorting the cropping limits in ascending order
xdat_lims       = sort(xdat_lims);
%% 1 - Cropping the data
% - x-axis indices
[~, xIndxL]     = min(abs(xdat(:) - xdat_lims(1)));
[~, xIndxU]     = min(abs(xdat(:) - xdat_lims(2)));
x_indx          = [xIndxL xIndxU];
% - cropping x-axis
xdat_crop       = xdat(x_indx(1):x_indx(2));
ydat_crop       = ydat(x_indx(1):x_indx(2));

end
