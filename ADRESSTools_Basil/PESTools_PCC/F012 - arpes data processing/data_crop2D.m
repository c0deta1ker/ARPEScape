function [xdat_crop, ydat_crop, ddat_crop] = data_crop2D(xdat, ydat, ddat, xdat_lims, ydat_lims)
% [xdat_crop, ydat_crop, ddat_crop] = data_crop2D(xdat, ydat, ddat, xdat_lims, ydat_lims)
%   This function crops 2D data along the x- and y-axis. The user can
%   define the limits along the x- and y-axis to be cropped.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%
%   IN:
%   -   xdat:           [1 x nX] or [nY x nX] array of the x-axis data
%   -   ydat:           [nY x 1] or [nY x nX] array of the y-axis data
%   -   ddat:           [nY x nX] array of the intensity data
% 	-   xdat_lims:      [1 x 2] row vector of x-axis limits
% 	-   ydat_lims:      [1 x 2] row vector of y-axis limits
%
%   OUT:
%   -   xdat_crop:  	[A x B] array of the cropped x-axis
%   -   ydat_crop:   	[A x B] array of the cropped y-axis
%   -   ddat_crop:  	[A x B] array of the cropped intensity data

%% Default parameters
maxLim = 1e4;
if nargin < 4; xdat_lims=[min(xdat(:)), max(xdat(:))]; ydat_lims=[min(ydat(:)), max(ydat(:))]; end
if nargin < 5; ydat_lims=[min(ydat(:)), max(ydat(:))]; end
if isempty(xdat_lims); xdat_lims=[min(xdat(:)), max(xdat(:))];  end
if isempty(ydat_lims); ydat_lims=[-1 1]*maxLim;  end
% - Sorting the cropping limits in ascending order
xdat_lims = sort(xdat_lims);
ydat_lims = sort(ydat_lims);
% - Initialising variables
xdat_crop = [];
ydat_crop = [];
ddat_crop = [];

%% 1 - Cropping the data
% - x-axis indices
if ~isempty(xdat)
    [~, xIndxL]     = min(abs(xdat(1,:) - xdat_lims(1)));
    [~, xIndxU]     = min(abs(xdat(1,:) - xdat_lims(2)));
    x_indx          = [xIndxL xIndxU];
end
% - y-axis indices
if ~isempty(ydat)
    [~, yIndxL]     = min(abs(ydat(:,1) - ydat_lims(1)));
    [~, yIndxU]     = min(abs(ydat(:,1) - ydat_lims(2)));
    y_indx          = [yIndxL yIndxU];
end
% - cropping x-axis
if ~isempty(xdat)
    if size(xdat, 1) > 1 && isempty(ydat); 	xdat_crop = xdat(:, x_indx(1):x_indx(2));
    elseif size(xdat, 1) > 1;               xdat_crop = xdat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
    else;                                   xdat_crop = xdat(1,x_indx(1):x_indx(2));
    end
end
% - cropping y-axis
if ~isempty(ydat)
    if size(ydat, 2) > 1 && isempty(xdat); 	ydat_crop = ydat(y_indx(1):y_indx(2), :);
    elseif size(ydat, 2) > 1;               ydat_crop = ydat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
    else;                                   ydat_crop = ydat(y_indx(1):y_indx(2),1);
    end
end
% - cropping the data
if isempty(ydat);       ddat_crop = ddat(:,x_indx(1):x_indx(2));
elseif isempty(xdat);   ddat_crop = ddat(y_indx(1):y_indx(2),:);
else;                   ddat_crop = ddat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
end
end
